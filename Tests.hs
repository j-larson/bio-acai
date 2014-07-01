module Main where

import Test.Hspec
import Bio.Sequence as S
import Data.Clustering.Hierarchical
import Global
import Output

--
-- test Output module
-- 

leaf :: String -> Int -> Dendrogram (AnnSeq a)
leaf fileName seqNumber = 
  let header = S.fromStr "Sequence header"
      content = S.fromStr "GATTACA"
      dummySequence = Seq header content Nothing
  in Leaf AnnSeq{seqFile=fileName, seqNum=seqNumber, seqData=dummySequence}

clusters :: [(Int, Dendrogram (AnnSeq a))]
clusters = [(1, Branch 0.8 (Branch 0.5 (leaf "A" 1) (leaf "B" 5)) (leaf "C" 8)),
            (2, Branch 0.8 (Branch 0.5 (leaf "C" 3) (leaf "A" 6)) (leaf "A" 7)),
            (3, leaf "B" 4),
            (4, leaf "A" 2)]

clustRes :: String
clustRes =  "cluster_1 A:1,B:5,C:8\n" ++
            "cluster_2 C:3,A:6,A:7\n" ++
            "cluster_3 B:4\n" ++
            "cluster_4 A:2\n"

matrixRes :: String
matrixRes =  ",A,B,C\n" ++
             "cluster_1 1,1,1\n" ++
             "cluster_2 1,0,1\n" ++
             "cluster_3 0,1,0\n" ++
             "cluster_4 1,0,0\n"

main :: IO ()
main = hspec $ do
  describe "Output" $ do
    it "outputs clusters" $ do
      clustersOutput clusters `shouldBe` clustRes
    it "output the incidence matrix" $ do
      matrixOutput clusters ["A", "B", "C"] `shouldBe` matrixRes
