#Zixuan Qiu s2777279
#Team member contribution:
#Zixuan Qiu did the pre-processing of vector 'a' and create vector 'b' which contains around 1000 most common words from cleaned 'a' (section 4 and 5). 
#setwd("C:/Users/cheryl/AllRepoProjects/Firstrepogroup")
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

#6a map the full token vector a to vocabulary b to produce M1
mlag <- 4
M1 <- match(a, b)

#6b
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

#7 build next.word function
next.word <- function(key, M, M1, w = rep(1, ncol(M)-1)) {
  key_len <- length(key)#the length of the key sequence
  
  #If the key is longer than mlag, keep the last mlag tokens
  if (key_len > mlag) {
    key <- tail(key, mlag)
    key_len <- mlag
  }
  #initialize  vectors for candidate next tokens and their probabilities
  u_all <- c()
  p_all <- c()
  
for (m in 1:key_len) {
  subM <- M[, m:key_len, drop = FALSE] #the submatrix that contains relevant columns
  
  ii <- colSums(!(t(subM) == key[m:key_len]))
  match_rows <- which(ii == 0 & is.finite(ii))
  
  if (length(match_rows) > 0) {
    u <- M[match_rows, key_len+1, drop = TRUE]
    u <- u[!is.na(u)]#remove missing values
    
    #remove punctuation tokens from candidates bcause M and M1 still base on original b with punctuations
    #words_u <- b[u]
    #valid_idx <- words_u %in% b_no_punct
    #u <- u[valid_idx]
    
    if (length(u) > 0) {
      #if candidate exists, assign corresponding probabilities proportional to the weight
      prob <- rep(w[m] / length(u), length(u))
      u_all <- c(u_all, u)
      p_all <- c(p_all, prob)
    }
  }
}
  if (length(u_all) == 0) {
    #if no candidates were found, go back to sampling from the global frequency
    tab <- table(M1)#count the word frequency
    valid_tokens <- as.integer(names(tab))[b[as.integer(names(tab))] %in% b_no_punct]
    probs <- as.numeric(tab)[b[as.integer(names(tab))] %in% b_no_punct]
    probs <- probs / sum(probs)
    return(sample(valid_tokens, 1, prob = probs))
  }
  
  
  return(sample(u_all, 1, prob = p_all))
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

#9 simulate sentences until full stop
simulate_sentence <- function(start_token, M, M1, b, mlag, w = rep(1, mlag)) {
  sentence_tokens <- as.integer(start_token)
  repeat {
    # only uses the last n=mlag tokens
    current_sequence <- tail(sentence_tokens, n = mlag)
    next_token <- next.word(current_sequence, M, M1, w)
    sentence_tokens <- c(sentence_tokens, next_token)
    #compare words
    if (b[sentence_tokens[length(sentence_tokens)]] == ".") break
  }
  return(sentence_tokens)
}
st <- select_start_token(M1,b)
tokens <- simulate_sentence(st, M, M1, b, mlag)
sentence_words <- b[tokens]
sentence <- paste(sentence_words, collapse = " ")
print(sentence)
