{-# LANGUAGE DeriveDataTypeable, RecordWildCards, StandaloneDeriving #-}

import Bio.Alignment.AlignData                  (toStrings)
import Bio.Alignment.AAlign             as A
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

deriving instance Data Linkage
deriving instance Typeable Linkage

data Args = Args { dataDir  :: String 
                 , minFrac  :: Double
                 , minScore :: Int
                 , linkage  :: Linkage
                 }
            deriving (Show, Data, Typeable)

opts :: Args
opts  = Args { dataDir  = def &= opt "." &= typ "DIR"
                 &= help "Directory to search for .faa files; default \".\""
             -- how to show the defaults automatically in usage message?
             , minFrac  = def &= opt (0.0 :: Double)
             , minScore = def &= opt (50  :: Int)      -- ???
             , linkage  = SingleLinkage &= opt SingleLinkage }
             &= summary "Sequence clustering based on local alignment scoring"

ident :: AnnSeq a -> String
ident s = seqFile s ++ ":" ++ show (seqNum s)

--readAmino :: FilePath -> IO [Sequence Amino]
--readAmino  = fmap (map castToAmino) . readFasta

infinity :: Double
infinity  = 1/0.0

gapPen :: (Int, Int)
gapPen  = (-5,-5)

trim :: FilePath -> FilePath
trim  = reverse . drop 4 . reverse

untrim :: FilePath -> FilePath
untrim  = (++ ".faa")

d :: Double -> Sequence a -> Sequence a -> Double
d minFrac x y
  | frac >= minFrac = 1.0 / fromIntegral score
  | otherwise = infinity
  where
    frac = fromIntegral tot / fromIntegral (length s1)
    tot  = sum (zipWith (\a b -> fromEnum (a == b && a /= '-')) s1 s2)
    (s1,s2)      = toStrings seqs
    (score,seqs) = A.local_align M.blosum62 gapPen x y

main :: IO ()
main  = withFile "clustering" WriteMode $ \cfile ->
        withFile "matrix"     WriteMode $ \mfile -> do
  Args {..}  <- cmdArgs opts
  fs <- map trim . filter (".faa" `isSuffixOf`)
          <$> getDirectoryContents dataDir
  seqss <- mapM (\f -> zipWith (AnnSeq f) [1..] <$> F.readFasta (untrim f)) fs
  let clusters = zip [1::Int ..] $
          C.dendrogram linkage (concat seqss) (d minFrac `on` seqData)
        `cutAt` (1.0 / fromIntegral minScore)

  forM_ clusters $ \(n,c) ->
    hPutStrLn cfile $
      "cluster_" ++ show n ++ " " ++ intercalate "," (map ident (elements c))

  hPutStrLn mfile $ "," ++ intercalate "," fs
  forM_ clusters $ \(n,c) -> do
    hPutStr   mfile $ "cluster_" ++ show n ++ " "
    hPutStrLn mfile $ intercalate "," $
      map (show . fromEnum . (`elem` map seqFile (elements c))) fs
