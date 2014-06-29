{-# LANGUAGE DeriveDataTypeable, RecordWildCards, StandaloneDeriving #-}

import Bio.Alignment.AlignData                  (toStrings)
import Bio.Alignment.SAlign             as A
import Bio.Alignment.Matrices           as M
import Bio.Sequence
import Bio.Sequence.Fasta               as F
import Control.Applicative
import Control.Monad
import Data.Clustering.Hierarchical     as C
import Data.Data
import Data.Function                            (on)
import Data.List                                (intercalate, isSuffixOf)
import System.Console.CmdArgs
import System.Directory                         (getDirectoryContents)
import System.IO

-- type of sequences annotated with source file information
data AnnSeq a = AnnSeq { seqFile :: FilePath
                       , seqNum  :: Int -- position in file
                       , seqData :: Sequence a }

-- data LinkageOpt = Single | Complete | Upga deriving (Show, Data, Typeable) -- rename

deriving instance Data Linkage
deriving instance Typeable Linkage

data Args = Args { dataDir  :: String 
                 , minFrac  :: Double
                 , minScore :: Int
                 , linkage  :: Linkage
                 -- , linkage  :: LinkageOpt
                 }
            deriving (Show, Data, Typeable)

-- TODO bring in Data.Default ??
opts :: Args
opts  = Args { dataDir  = "." &= opt "."  -- ??? redundant; how to show the defaults automatically in help?
                              &= help "Directory to search for .faa files; default ."
             , minFrac  = 0.0 &= opt (0.0 :: Double) -- ???
             , minScore = 50  &= opt (50  :: Int)    -- ???
             , linkage  = SingleLinkage &= opt SingleLinkage }     -- ???
             &= summary "Sequence clustering based on local alignment scoring"

ident :: AnnSeq a -> String
ident s = seqFile s ++ ":" ++ show (seqNum s)

--readAmino :: FilePath -> IO [Sequence Amino]
--readAmino  = fmap (map castToAmino) . readFasta

infinity :: Double
infinity  = 1/0.0

gapPen :: Int
gapPen  = -5

trim :: FilePath -> FilePath
trim  = reverse . drop 4 . reverse

untrim :: FilePath -> FilePath
untrim  = (++ ".faa")

main :: IO ()
main  = withFile "clustering" WriteMode $ \cfile ->
        withFile "matrix"     WriteMode $ \mfile -> do
  Args {..}  <- cmdArgs opts
  fs <- map trim . filter (".faa" `isSuffixOf`) <$> getDirectoryContents dataDir
  seqss <- mapM (\f -> zipWith (AnnSeq f) [1..] `fmap` F.readFasta (untrim f)) fs
  -- let linkage' = case linkage of
  --       Single   -> C.SingleLinkage
  --       Complete -> C.CompleteLinkage
  --       Upga     -> C.UPGMA
  let d x y
        | frac >= minFrac =
        -- this recomputation of the score is spurious as we know the alignment
        -- so we can find a nicer way of extracting it:
            1.0 / fromIntegral (A.local_score M.blosum62 gapPen x y)
        | otherwise = infinity
        where
          frac    = fromIntegral tot / fromIntegral (length s1)
          tot     = sum (zipWith (\a b -> fromEnum (a == b && a /= '-')) s1 s2)
          (s1,s2) = toStrings $ A.local_align M.blosum62 gapPen x y
      clusters = zip [1::Int ..] $
          C.dendrogram linkage (concat seqss) (d `on` seqData)
        `cutAt` (1.0 / fromIntegral minScore)

  forM_ clusters $ \(n,c) ->
    hPutStrLn cfile $
      "cluster_" ++ show n ++ " " ++ intercalate "," (map ident (elements c))

  hPutStrLn mfile $ "," ++ intercalate "," fs
  forM_ clusters $ \(n,c) -> do
    hPutStr   mfile $ "cluster_" ++ show n ++ " "
    hPutStrLn mfile $ intercalate "," $
      map (show . fromEnum . (`elem` map seqFile (elements c))) fs
