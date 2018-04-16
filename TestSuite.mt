(* Wolfram Language Test file *)

dir = "Test";
files = {};
TestSuite[FileNameJoin[{dir, #}] & /@ files];
