(ns pubmed-search.core
  (:gen-class) ;; This allows compiling as a binary
  (:require hashp.core
            [clojure.data.json :as json]
            [clojure.data.xml :as xml]
            [clojure.edn :as edn]
            [clojure.java.io :as io]
            [clojure.string :as str]
            [donut.system :as ds]
            [hato.client :as hc]
            [honey.sql :as sql]
            [json-parse.core :as json-parse]
            [next.jdbc :as jdbc]
            [next.jdbc.result-set :as result-set]
            [org.httpkit.server :as server]
            [reitit.core :as re]
            [reitit.ring :as rr]
            [reitit.ring.middleware.parameters :refer [parameters-middleware]]
            [salmon.signal :as sig]))

(defonce state (atom nil))

(def jdbc-opts
  {:builder-fn result-set/as-kebab-maps})

(defn execute! [sqlite sqlmap]
  (jdbc/execute! (:datasource sqlite sqlite) (sql/format sqlmap) jdbc-opts))

(defn execute-one! [sqlite sqlmap]
  (jdbc/execute-one! (:datasource sqlite sqlite) (sql/format sqlmap) jdbc-opts))

;; Rate limit is 10 requests per second per API key.
;; We try to limit requests to once per 100 ms to stay under on average.
;; If we still get a 429, hold off on requests for 1000 ms and retry once.
(defonce last-request (atom (System/nanoTime)))

(defn api-get [url {:keys [throw-exceptions?] :as opts}]
  (locking last-request
    (let [diff-ms (quot (- (System/nanoTime) @last-request) 1000000)]
      (when (> 100 diff-ms)
        (Thread/sleep diff-ms))
      (reset! last-request (System/nanoTime))))
  (let [{:keys [status] :as response} (hc/get url (assoc opts :throw-exceptions? false))]
    (cond
      (= 429 status)
      (do
        (locking last-request
          (Thread/sleep 1000)
          (reset! last-request (System/nanoTime)))
        (hc/get url opts))

      (<= 200 status 300)
      response

      (false? throw-exceptions?)
      (throw (ex-info (str "Unexpected status: " status)
                      {:opts opts
                       :response response
                       :url url}))

      :else response)))

(defn tag-val [tag-name content]
  (some #(when (= tag-name (:tag %)) %) content))

(defn get-search-results-api [{:keys [api-key hato] :as opts} query & [retstart]]
  (let [docs (-> (api-get
                  "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
                  {:http-client (:client hato)
                   :query-params
                   {:db "pubmed"
                    :key api-key
                    :retmax 100
                    :retstart (or retstart 0)
                    :term query}
                   :timeout 10000})
                 :body
                 xml/parse-str
                 :content)
        ct (-> (tag-val :Count docs) :content first parse-long)
        ret-start (-> (tag-val :RetStart docs) :content first parse-long)
        id-list (-> (tag-val :IdList docs) :content)
        next-start (+ ret-start (count id-list))]
    (concat
     (mapcat :content id-list)
     (when (and (seq docs) (> ct next-start))
       (get-search-results-api opts query next-start)))))

(defn get-search-results-cache
  "Returns cached search results. Returns nil if the query is not cached.
   Returns () if the query is cached but returned no results."
  [{:keys [cache-ttl-sec sqlite]} query]
  (let [results (->> {:select [:pmid]
                      :from :pubmed-search-result
                      :where [:and
                              [:= :term query]
                              [:= :created
                               {:select [:%max.created]
                                :from :pubmed-search-result
                                :where [:and
                                        [:= :term query]
                                        [:>= :created
                                         [:datetime "now" (str "-" cache-ttl-sec " seconds")]]]}]]}
                     (execute! sqlite)
                     (map :pubmed-search-result/pmid))]
    (if (= [nil] results)
      ()
      (seq results))))

(defn put-search-results-cache [{:keys [sqlite]} query results]
  (jdbc/with-transaction [tx (:datasource sqlite)]
    (execute-one!
     tx
     {:delete-from :pubmed-search-result
      :where [:= :term query]})
    (execute-one!
     tx
     {:insert-into :pubmed-search-result
      :values
      (if (seq results)
        (map #(do {:pmid % :term query}) results)
        [{:pmid nil :term query}])}))
  results)

(defn get-search-results [opts query]
  (or
   (get-search-results-cache opts query)
   (->> (get-search-results-api opts query)
        (put-search-results-cache opts query))))

(defn get-fetches-api [{:keys [api-key hato]} pmids]
  (-> (api-get
       "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
       {:http-client (:client hato)
        :query-params
        {:db "pubmed"
         :id (str/join "," pmids)
         :key api-key
         :retmode "text"
         :rettype "xml"}
        :timeout 10000})
      :body
      xml/parse-str
      :content
      (->> (filter #(= :PubmedArticle (:tag %))))))

(defn get-fetches-cache
  "Returns cached PubMed fetches."
  [{:keys [cache-ttl-sec sqlite]} pmids]
  (->> {:select [:fetch]
        :from :pubmed-fetch
        :where [:and
                [:in :pmid pmids]
                [:>= :created
                 [:datetime "now" (str "-" cache-ttl-sec " seconds")]]]}
       (execute! sqlite)
       (map (comp xml/parse-str :pubmed-fetch/fetch))))

(defn fetch-id [fetch]
  (->> fetch
       :content
       (tag-val :MedlineCitation)
       :content
       (tag-val :PMID)
       :content
       first))

(defn put-fetches-cache [{:keys [sqlite]} fetches]
  (when (seq fetches)
    (jdbc/with-transaction [tx (:datasource sqlite)]
      (execute-one!
       tx
       {:delete-from :pubmed-fetch
        :where [:in :pmid (map fetch-id fetches)]})
      (execute-one!
       tx
       {:insert-into :pubmed-fetch
        :values
        (map
         #(do {:pmid (fetch-id %) :fetch (xml/emit-str %)})
         fetches)})))
  fetches)

(defn get-fetches-seq [{:keys [chunk-size] :as opts} cached pmid-set]
  (lazy-seq
   (if (seq cached)
     (let [[head & more] cached]
       (cons head (get-fetches-seq opts more (disj pmid-set (fetch-id head)))))
     (when (seq pmid-set)
       (let [chunk (take chunk-size pmid-set)]
         (concat
          (->> (get-fetches-api opts chunk)
               (put-fetches-cache opts))
          (get-fetches-seq opts nil (reduce disj pmid-set chunk))))))))

(defn get-fetches [opts pmids]
  (when (seq pmids)
    (get-fetches-seq opts (get-fetches-cache opts pmids) (set pmids))))

(defn get-config [filename]
  (if-let [resource (if (.exists (io/file filename))
                      (io/file filename)
                      (io/resource filename))]
    (with-open [reader (-> resource io/reader java.io.PushbackReader.)]
      (try
        (edn/read reader)
        (catch Exception e
          (throw
           (ex-info (str "Error parsing EDN in config file \"" filename
                         \" ": " (.getMessage e))
                    {:filename filename}
                    e)))))
    (throw
     (ex-info (str "Config file not found: \"" filename "\"")
              {:filename filename}))))

(defn config-component []
  #::ds{:start (fn [{::ds/keys [config]}]
                 (get-config (:filename config)))
        :stop (constantly nil)
        :config {:filename (ds/local-ref [:env :config-file])}})

(defn hato-component [config]
  {::ds/config config
   ::ds/start
   (fn [{::ds/keys [config instance]}]
     (let [{:keys [client]} instance]
       (if client
         instance
         (assoc instance :client
                (hc/build-http-client config)))))
   ::ds/stop
   (fn [{::ds/keys [instance]}]
     (if (:client instance)
       (assoc instance :client nil)
       instance))})

(def sql-create-schema
  [; A null pmid is used to indicate empty search results
   "CREATE TABLE IF NOT EXISTS pubmed_search_result (
      id INTEGER PRIMARY KEY,
      term TEXT NOT NULL,
      created TEXT NOT NULL DEFAULT CURRENT_TIMESTAMP,
      pmid TEXT
    ) STRICT"
   "CREATE INDEX IF NOT EXISTS term_index ON pubmed_search_result (term)"
   "CREATE INDEX IF NOT EXISTS created_index ON pubmed_search_result (created)"
   "CREATE INDEX IF NOT EXISTS pmid_index ON pubmed_search_result (pmid)"
   "CREATE TABLE IF NOT EXISTS pubmed_fetch (
      id INTEGER PRIMARY KEY,
      pmid TEXT NOT NULL,
      created TEXT NOT NULL DEFAULT CURRENT_TIMESTAMP,
      fetch TEXT NOT NULL
    ) STRICT"
   "CREATE INDEX IF NOT EXISTS pmid_index ON pubmed_fetch (pmid)"
   "CREATE INDEX IF NOT EXISTS created_index ON pubmed_fetch (created)"])

(defn sqlite-component [config]
  {::ds/config config
   ::ds/start
   (fn [{::ds/keys [config instance]}]
     (let [{:keys [datasource]} instance
           db-spec (select-keys config [:dbname :dbtype])]
       (if datasource
         instance
         (let [datasource (jdbc/get-datasource db-spec)]
           (doseq [s sql-create-schema]
             (execute! datasource {:raw s}))
           (assoc instance :datasource datasource)))))
   ::ds/stop
   (fn [{::ds/keys [instance]}]
     (if (:client instance)
       (assoc instance :client nil)
       instance))})

(defn article-uri [pmid]
  (str "https://www.ncbi.nlm.nih.gov/pubmed/" pmid))

(defn doc-abstract [doc-json]
  (-> doc-json
      :PubmedArticle
      :MedlineCitation
      :Article
      (as->
       $
       (cond
         (map? $) (:Abstract $)
         (seq? $) (map :Abstract $))
        (cond
          (map? $) (:AbstractText $)
          (seq? $) (map :AbstractText $))
        (cond
          (string? $) [$]
          (seq? $) (keep identity $))
        (when $
          (str/join "\n\n" $)))))

(defn doc-title [doc-json]
  (-> doc-json
      :PubmedArticle
      :MedlineCitation
      :Article
      (as->
       $
       (cond
         (map? $) (:ArticleTitle $)
         (seq? $) (map :ArticleTitle $))
        (cond
          (string? $) [$]
          (seq? $) (keep identity $))
        (when $
          (str/join "\n\n" $)))))

(defn pubmed->srvc [doc-xml]
  (let [doc-json (json-parse/xml->json [doc-xml])]
    {:data (-> doc-json
               (assoc :abstract (doc-abstract doc-json)
                      :title (doc-title doc-json)))
     :uri (article-uri (fetch-id doc-xml))
     :type "document"}))

(defn get-home [{:keys [query-params]} opts]
  (let [q (get query-params "q")]
    (if (seq q)
      (->> (get-search-results opts q)
           (get-fetches opts)
           (map pubmed->srvc)
           (map json/write-str)
           (interpose "\n")
           (hash-map :status 200 :body))
      {:status 200
       :body ""})))

(defn routes [opts]
  [["/" {:get #(get-home % opts)
         :middleware [parameters-middleware]}]])

(defn http-server-component [config]
  #::ds{:config config
        :start (fn [{:keys [::ds/config]}]
                 (let [server (server/run-server
                               (-> (routes config)
                                   rr/router
                                   rr/ring-handler)
                               {:legacy-return-value? false
                                :port (:port config)})]
                   {:port (server/server-port server)
                    :server server}))
        :stop (fn [{::ds/keys [instance]}]
                @(server/server-stop! (:server instance))
                nil)})

(defn system [env]
  {::ds/defs
   {:search
    {:config (config-component)
     :env {::ds/start (constantly env)}
     :hato (hato-component (ds/local-ref [:config :hato]))
     :http-server (http-server-component
                   {:api-key (ds/local-ref [:env :api-key])
                    :cache-ttl-sec (ds/local-ref [:config :cache-ttl-sec])
                    :chunk-size (ds/local-ref [:config :chunk-size])
                    :hato (ds/local-ref [:hato])
                    :port (ds/local-ref [:config :port])
                    :sqlite (ds/local-ref [:sqlite])})
     :sqlite (sqlite-component {:dbname (ds/local-ref [:config :db])
                                :dbtype "sqlite"})}}})

;; Not thread-safe. For use by -main and at REPL
(defn start! []
  (let [env {:api-key (System/getenv "EUTILS_API_KEY")
             :config-file "pubmed-search-config.edn"
    (swap! state #(sig/signal! (or % (system env)) ::ds/start))))

;; Not thread-safe. For use by -main and at REPL
(defn stop! []
  (swap! state #(when % (sig/signal! % ::ds/stop) nil)))

(defn -main []
  (start!)
  (Thread/sleep Long/MAX_VALUE))

(comment
  (do (stop!) (start!) nil))
