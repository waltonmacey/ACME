#!/bin/bash


branch=$1
strD="These files were deleted"
strA="These files were added"
str1="This is the original assignment variable"
str2="This is the new assignment variable"

#git diff --numstat $branch | while read add del file
#do echo $file
#done

printf "\n"
echo $str1
git diff $branch | grep -m 4 '^-    ' | awk -F'=>|%' '{print $2;}' | awk '!seen[$0]++'

printf "\n"
echo $str2
git diff $branch | grep -m 4 '^+    ' | awk -F'=>|%' '{print $2;}' | awk '!seen[$0]++'


#git diff $branch | grep -m 4 '^-    ' | sed 's/.*=> \(.*\)\%/\1/'

printf "\n"
git diff --shortstat $branch
printf "\n"

echo $strA
git diff --name-status $branch | grep '^A' 
echo $strD
git diff --name-status $branch | grep '^D' 


