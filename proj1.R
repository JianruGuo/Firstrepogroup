library(debug)
mtrace(split_punct)
setwd("C:/Users/cheryl/AllRepoProjects/Firstrepogroup")
a <- scan("shakespeare.txt",what="character",skip=83,nlines=196043-83,
          fileEncoding="UTF-8")

#4a remove stage directions
#use grep function to locate words that contain "[" in vector a
open_bracket_idx <- grep("\\[",a)

#create a logical vector to store the property of all words for later removing
remove_idx <- rep(FALSE,length(a))

#loop through the index of open brackets  
for (i in open_bracket_idx) {
  max_look_ahead <- min(100, length(a) - i) #make sure not exceed the text limit
  look_ahead <- i + (1:max_look_ahead)
  close_bracket_idx <- grep("\\]",a[look_ahead]) #use grep function to locate all "]" within next 100 words, return index of the sub-vector
  if (length(close_bracket_idx) > 0) { #make sure there exist close brackets within next 100 words
    first_close_idx <- close_bracket_idx[1]  #find the first corresponding ']' in the sub-vector that matches '['
    actual_end_idx <- i + first_close_idx #the actual index in a
    remove_idx[i:actual_end_idx] <- TRUE #build a record of all words corresponding to stage directions
  }
}
#take the opposite logical value of all words, keep the words which are TRUE, i.e. removing stage directions
a_clean <- a[!remove_idx]


#4b remove fully upper case words and numerals
words <- a_clean
#create logical vector; TRUE for fully upper case words except "I" and "A", FALSE for others
is_upper <- (words == toupper(words)) & !(words %in% c("I", "A"))
#use grepl function to return a logical vector, TRUE for numerals; "^[0-9]+$" is the regular expression of numbers
is_number <- grepl("^[0-9]+$", words)
#remove fully upper case words and numerals simultaneously to avoid unmatched index when remove them separately
a_remain <- words[!(is_upper | is_number)] 


#4c
#use gsub function to remove "-" and "_", i.e. replace "-" and "_" with ""
a <- gsub("[-_]","",a_remain)


#4d separate punctuation as new entries
split_punct <- function(a, punct_vec) {
  for (p in punct_vec){
    punct_loc <- grep(p, a, fixed = TRUE) #index of words containing punctuation
    if (length(punct_loc) > 0) { #make sure punctuation exists
      a_clean <- gsub(p,"",a,fixed = TRUE) #remove punctuation from words
      new_length <- length(a) + length(punct_loc) 
      result <- character(new_length) #create new vector to store the result
      punct_insert_positions <- punct_loc + 1:length(punct_loc) #location where we insert punctuation
      result[punct_insert_positions] <- p #insert all punctuation
      result[-punct_insert_positions] <- a_clean #insert all other words
      a <- result
    }
  }
  return(a)
}

#4e
punctuations <- c(",", ".", ";", "!", ":", "?")
a <- split_punct(a,punctuations) #Assign the result back to "a"

#4f convert all words to lower case
a <- tolower(a)


#5a use unique function to create vector b of unique words from "a"
b_unique <- unique(a)

#5b use match function to find the index of unique words in "a"
unique_idx <- match(a,b_unique)

#5c use tabulate function to count the number of unique words
word_counts <- tabulate(unique_idx, nbins = length(b_unique))

#5d use rank function to rank words from the most common to least common
ranks <- rank(-word_counts,ties.method = "min")
#take about 1000 most common words and save as vector "b"
b <- b_unique[ranks <= 1000]

#6a
mlag <- 4

#6b map the full token vector a to vocabulary b to produce M1
M1 <- match(a, b)

#6c
#dimension of matrix锛?(n - mlag) x (mlag + 1)
#Define build_M, which takes M1 and mlag and returns the shifted (lag) matrix M锛?
build_M <- function(M1, mlag) {
  n <- length(M1)
  #initialize with NA_integer_ (integer NA), with n-mlag rows and mlag+1 columns
  out <- matrix(NA_integer_, n - mlag, mlag + 1)
  for (j in 0:mlag) {
    out[, j + 1] <- M1[(1 + j):(n - mlag + j)]
  }
  out
}
M <- build_M(M1, mlag)

