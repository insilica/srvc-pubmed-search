{
  description = "pubmed-search";

  inputs = {
    nixpkgs.url = "github:nixos/nixpkgs/nixos-22.05";
    nixpkgs-unstable.url = "github:nixos/nixpkgs/nixpkgs-unstable";
    flake-utils.url = "github:numtide/flake-utils";
    clj-nix = {
      url = "github:jlesquembre/clj-nix";
      inputs.flake-utils.follows = "flake-utils";
      inputs.nixpkgs.follows = "nixpkgs";
    };
    srvc.url = "github:insilica/rs-srvc";
  };
  outputs =
    { self, nixpkgs, nixpkgs-unstable, flake-utils, clj-nix, srvc, ... }@inputs:
    flake-utils.lib.eachDefaultSystem (system:
      with import nixpkgs { inherit system; };
      let
        pkgs-unstable = import nixpkgs-unstable { inherit system; };
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
          name = "pubmed-http-server";
          main-ns = "pubmed-search.core";
          jdkRunner = pkgs.jdk17_headless;
        };
        # Strips out unused JDK code for a smaller binary
        http-server = cljpkgs.customJdk { cljDrv = http-server-bin; };
        http-server-image = dockerTools.buildLayeredImage {
          name = "insilica/pubmed-http-server";
          tag = "latest";
          config = {
            Cmd = clj-nix.lib.mkCljCli {
              jdkDrv = http-server;
              java-opts = [ "-Dclojure.compiler.direct-linking=true" ];
            };
            ExposedPorts = { "7881" = { }; };
          };
          contents = [ sqlite ];
          extraCommands = "mkdir -m 0777 tmp";
        };
      in {
        packages = {
          inherit http-server http-server-image pubmed-search;
          default = pubmed-search;
        };
        devShells.default = mkShell {
          buildInputs = [
            clj-nix.packages.${system}.deps-lock
            clojure
            pkgs-unstable.flyctl
            my-python
            pubmed-search
            srvc.packages.${system}.default
            yapf
          ];
        };
      });
}
