{
  description = "pubmed-search";

  inputs = {
    nixpkgs.url = "github:nixos/nixpkgs/nixos-22.05";
    flake-utils.url = "github:numtide/flake-utils";
    srvc.url = "github:insilica/rs-srvc";
  };
  outputs = { self, nixpkgs, flake-utils, srvc, ... }@inputs:
    flake-utils.lib.eachDefaultSystem (system:
      with import nixpkgs { inherit system; };
      let
        my-python = python3.withPackages (ps: with ps; [ ps.biopython ]);
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
      in {
        packages = {
          inherit pubmed-search;
          default = pubmed-search;
        };
        devShells.default = mkShell {
          buildInputs =
            [ my-python pubmed-search srvc.packages.${system}.default yapf ];
        };
      });
}
