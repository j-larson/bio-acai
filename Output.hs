module Output
( clustersOutput
, matrixOutput
) where

import Data.Clustering.Hierarchical
import Data.List
import Global

-- Given a numbered set of clusters, produce an output string like this:
--
-- cluster_1 NC_012654:1,NC_012655:1
-- cluster_2 NC_012654:2,NC_012655:2
-- cluster_3 NC_012655:3
-- cluster_4 NC_012654:3
--
clustersOutput :: [(Int, Dendrogram (AnnSeq a))] -> String
clustersOutput clusters = unlines $ map cline clusters
  
cline :: (Int, Dendrogram (AnnSeq a)) -> String
cline (n,c) = 
  "cluster_" ++ show n ++ " " ++ intercalate "," (map ident (elements c))

ident :: AnnSeq a -> String
ident s = seqFile s ++ ":" ++ show (seqNum s)

-- Given a numbered set of clusters and a sequence of input files,
-- produce an output string like this:
-- ,NC_012654,NC_012655
-- cluster_1,1,1
-- cluster_2,1,1
-- cluster_3,0,1
-- cluster_4,1,0
matrixOutput :: [(Int, Dendrogram (AnnSeq a))] -> [String] -> String
matrixOutput clusters fs = 
  let header = "," ++ intercalate "," fs
      mains = map (matrixLine fs) clusters
  in unlines (header : mains)

matrixLine :: [String] -> (Int, Dendrogram (AnnSeq a)) -> String
matrixLine fs (n,c) =
  let intro = "cluster_" ++ show n ++ " " 
      digits = map (show . fromEnum . (`elem` map seqFile (elements c))) fs
  in intro ++ intercalate "," digits
