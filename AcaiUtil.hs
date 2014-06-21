module AcaiUtil
( upperCase
, lowerCase
) where

import Data.Char

upperCase :: String -> String
upperCase s = map toUpper s 

lowerCase :: String -> String
lowerCase s = map toLower s 
