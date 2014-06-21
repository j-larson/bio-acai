module Main where

import AcaiUtil
import Test.Hspec

main :: IO ()
main = hspec $ do
  describe "AcaiUtils" $ do
    it "converts to upper case" $ do
      upperCase "simple string" `shouldBe` "SIMPLE STRING"

    it "converts to lower case" $ do
      lowerCase "STUFF" `shouldBe` "stuffa" 
