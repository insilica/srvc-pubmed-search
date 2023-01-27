{
  description = "pubmed-search";

  inputs = {
    nixpkgs.url = "github:nixos/nixpkgs/nixos-22.05";
    flake-utils.url = "github:numtide/flake-utils";
    clj-nix = {
      url = "github:jlesquembre/clj-nix";
      inputs.flake-utils.follows = "flake-utils";
      inputs.nixpkgs.follows = "nixpkgs";
    };
    srvc.url = "github:insilica/rs-srvc";
  };
  outputs = { self, nixpkgs, flake-utils, clj-nix, srvc, ... }@inputs:
    flake-utils.lib.eachDefaultSystem (system:
      with import nixpkgs { inherit system; };
      let
        cljpkgs = clj-nix.packages."${system}";
        my-python = python3.withPackages (ps: with ps; [ biopython xmltodict ]);
        pubmed-search = stdenv.mkDerivation {
          pname = "pubmed-search";
          version = "0.1.0";
          src = ./src;
          buildInputs = [ my-python ];
          installPhase = ''
            mkdir -p $out/bin
            cp pubmed-search.py $out/bin/pubmed-search
          '';
        };
        http-server-bin = cljpkgs.mkCljBin {
          projectSrc = ./.;
          name = "srvc-pubmed-search-http-server";
          main-ns = "pubmed-search.core";
          jdkRunner = pkgs.jdk17_headless;
        };
        # Strips out unused JDK code for a smaller binary
        http-server = cljpkgs.customJdk {
          cljDrv = http-server-bin;
        };
      in {
        packages = {
          inherit http-server pubmed-search;
          default = pubmed-search;
        };
        devShells.default = mkShell {
          buildInputs = [
            clj-nix.packages.${system}.deps-lock
            clojure
            my-python
            pubmed-search
            srvc.packages.${system}.default
            yapf
          ];
        };
      });
}
