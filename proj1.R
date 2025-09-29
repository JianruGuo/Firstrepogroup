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

#7
next.word <- function(key, M, M1, w = rep(1, ncol(M) - 1)) {
  
  if (length(key) > mlag) {
    key <- key[(length(key) - mlag + 1):length(key)]# If the length of key is too long, use the last mlag words.
  }
  
  current_length <- length(key)
 
  all_candidates <- integer(0)# initialize candidate words and weights
  all_weights <- numeric(0)
  

  for (l in current_length:1) {
    mc <- mlag - l + 1
    me <- mlag  # match the end column
    subM <- M[, mc:me, drop = FALSE]# get submatrices of M
    
    #transpose M to make the row of M match the row of key
    key_subset <- key[(current_length - l + 1):current_length]
    comparison_matrix <- t(subM) == key_subset
    ii <- colSums(!comparison_matrix)
    
    matching_rows <- which(ii == 0 & is.finite(ii))#find rows that match key
    
    if (length(matching_rows) > 0) {
      #get next word token
      candidates <- M[matching_rows, mlag + 1]
      candidates <- candidates[!is.na(candidates)]
      
      if (length(candidates) > 0) {
        #calculate mixture weight
        weight <- w[mlag - l + 1] / length(candidates)
        weights <- rep(weight, length(candidates))
        all_candidates <- c(all_candidates, candidates)
        all_weights <- c(all_weights, weights)
      }
    }
  }
  
  
  if (length(all_candidates) == 0) {
    # if no candidate word is found
    valid_tokens <- M1[!is.na(M1)]
    if (length(valid_tokens) > 0) {
      return(sample(valid_tokens, 1))
    } else {
      return(1)  
    }
  }
  
  # using the sample function to sample one token from this distribution
  selected_index <- sample(1:length(all_candidates), 1, prob = all_weights)
  return(all_candidates[selected_index])
}


#8 select a single word token(but not punctuation) at random 
select_start_token <- function(M1, b, start_word = NULL) {
    # define punctuation marks in order to exclude them
    punct_marks <- c(",", ".", ";", "!", ":", "?")
    
    # get tokens that are non-NA and not punctuation
    valid_tokens <- M1[!is.na(M1) & !(b[M1] %in% punct_marks)]
    
    if (length(valid_tokens) == 0) {
      stop("No valid starting tokens found")
    }
    
    start_token <- sample(valid_tokens, 1)
    cat("Using random starting word: '", b[start_token], "' (token:", start_token, ")\n", sep = "")
    
    return(start_token)
}

#9 - Simulate sentences until full stop
simulate_sentence <- function(start_token, M, M1, b, mlag, w = rep(1, mlag)) {
  current_sequence <- start_token
  sentence_tokens <- c(start_token)
  
  # Continue until we reach a full stop or max iterations
  while (iteration < max_iterations) {
    iteration <- iteration + 1
    
    # get next word using the next_word function
    next_token <- next.word(current_sequence, M, M1, w)
    
    # add to our sentence
    sentence_tokens <- c(sentence_tokens, next_token)
    
    # update current sequence (keep only last mlag tokens)
    if (length(sentence_tokens) >= mlag) {
      current_sequence <- sentence_tokens[(length(sentence_tokens) - mlag + 1):length(sentence_tokens)]
    } else {
      current_sequence <- sentence_tokens
    }
    
    # check if we reached a full stop
    if (b[next_token] == ".") {
      cat("Full stop reached after", iteration, "iterations\n")
      break
    }
  }
    
    # Convert tokens back to words
    sentence_words <- b[sentence_tokens]
    
    # Format the sentence nicely
    formatted_sentence <- format_sentence(sentence_words)
    
    return(list(
      tokens = sentence_tokens,
      words = sentence_words,
      sentence = formatted_sentence,
      length = length(sentence_tokens)
    ))
}
  
  # main function to run the complete simulation
run_shakespeare_simulator <- function(M, M1, b, mlag = 4, start_word = NULL) {
  cat("=== Shakespeare Sentence Simulator ===\n")
  start_token <- select_start_token(M1, b, start_word)# select starting token
  result <- simulate_sentence(start_token, M, M1, b, mlag)# simulate sentence
    
  # print results
  cat("\n=== Generated Sentence ===\n")
  cat(result$sentence, "\n")
  cat("\nSentence length:", result$length, "words\n")
  cat("Word tokens:", paste(result$tokens, collapse = " "), "\n")
    
  return(result)
}

