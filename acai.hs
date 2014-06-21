import Data.Char (toUpper)

upperCase :: String -> String
upperCase s = map toUpper s 

main :: IO ()
main = do
  putStrLn "Hello world!"
  putStrLn $ upperCase "Welcome to Acai."
