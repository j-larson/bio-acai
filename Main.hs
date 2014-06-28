import Bio.Alignment.AlignData                  (toStrings)
import Bio.Alignment.SAlign             as A
import Bio.Alignment.Matrices           as M
import Bio.Sequence
import Bio.Sequence.Fasta               as F
import Control.Monad
import Data.Clustering.Hierarchical     as C
import Data.Function                            (on)
import Data.List                                (intercalate)
import System.Environment                       (getArgs)
import System.IO

-- type of sequences annotated with source file information
data AnnSeq a = AnnSeq { seqFile :: FilePath
                       , seqNum  :: Int -- position in file
                       , seqData :: Sequence a } deriving Show

ident :: AnnSeq a -> String
ident s = seqFile s ++ ":" ++ show (seqNum s)

--readAmino :: FilePath -> IO [Sequence Amino]
--readAmino  = fmap (map castToAmino) . readFasta

infinity :: Double
infinity  = 1/0.0

gapPen, minScore :: Int
gapPen   = -5
minScore = 50

minFrac :: Double
minFrac  = 0.8

main :: IO ()
main  = withFile "clustering" WriteMode $ \cfile ->
        withFile "matrix"     WriteMode $ \mfile -> do
  fs    <- getArgs
  seqss <- mapM (\f -> zipWith (AnnSeq f) [1..] `fmap` F.readFasta f) fs
  let d x y
        | frac >= minFrac =
        -- this recomputation of the score is spurious as we know the alignment
        -- so we can find a nicer way of extracting it:
            1.0 / fromIntegral (A.local_score M.blosum62 gapPen x y)
        | otherwise = infinity
        where
          frac    = fromIntegral tot / fromIntegral l
          tot     = sum (zipWith (\a b -> fromEnum (a == b && a /= '-')) s1 s2)
          (s1,s2) = toStrings $ A.local_align M.blosum62 gapPen x y
          l       = length s1
      clusters = zip [1::Int ..] $
          C.dendrogram SingleLinkage (concat seqss) (d `on` seqData)
        `cutAt` (1.0 / fromIntegral minScore)

  forM_ clusters $ \(n,c) ->
    hPutStrLn cfile $
      "cluster_" ++ show n ++ " " ++ intercalate "," (map ident (elements c))

  hPutStrLn mfile $ "," ++ intercalate "," fs
  forM_ clusters $ \(n,c) -> do
    hPutStr   mfile $ "cluster_" ++ show n ++ " "
    hPutStrLn mfile $ intercalate "," $
      map (show . fromEnum . (`elem` map seqFile (elements c))) fs
