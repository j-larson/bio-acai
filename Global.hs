module Global where

import Bio.Sequence (Sequence)

data AnnSeq a = AnnSeq { seqFile :: FilePath
                       , seqNum  :: Int -- position in file
                       , seqData :: Sequence a }
