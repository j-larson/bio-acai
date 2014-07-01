bio-acai
========

This is a tiny script for hierarchical clustering based on
local alignment scores.

To install, you need `ghc` and `cabal` installed.  
`cd` into the source directory and do

   cabal sandbox init   # optional
   cabal install --only-dependencies --allow-newer
   cabal build

This will create an executable at `./dist/build/acai/acai`.
You can get usage information with the `--help` flag.
