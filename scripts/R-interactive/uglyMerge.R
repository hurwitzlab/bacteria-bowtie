#uglyMerge####
#from http://stackoverflow.com/questions/10617377/merge-data-with-partial-match-in-r
uglyMerge <- function(df1, df2) {

  ## lower all strings to allow case-insensitive comparison
  lowerNames1 <- tolower(df1[, 1]);
  lowerNames2 <- tolower(df2[, 1]);

  ## split strings into single characters
  names1 <- strsplit(lowerNames1, "");
  names2 <- strsplit(lowerNames2, "");

  ## create the final dataframe
  mergedDf <- data.frame(name1=as.character(df1[,1]), name2=NA, number1=df1[,2], number2=NA, matches=0,stringsAsFactors=FALSE);

  ## store names of dataframe2 (to remember which strings have no match)
  toMerge <- df2[, 1];

  for (i in seq(along=names1)) {
    for (j in seq(along=names2)) {
      ## set minimal match to 4 or to string length
      minMatch <- min(4, length(names2[[j]]));

      ## find single matches
      matches <- names1[[i]] %in% names2[[j]];

      ## look for consecutive matches
      r <- rle(matches);

      ## any matches found?
      if (any(r$values)) {
        ## find max consecutive match
        possibleMatch <- r$value == TRUE;
        maxPos <- which(which.max(r$length[possibleMatch]) & possibleMatch)[1];

        ## store max conscutive match length
        maxMatch <- r$length[maxPos];

        ## to remove FALSE-POSITIVES (e.g. CCC and kcn) find
        ## largest substring
        start <- sum(r$length[0:(maxPos-1)]) + 1;
        stop <- start + r$length[maxPos] - 1;
        maxSubStr <- substr(lowerNames1[i], start, stop);

        ## all matching criteria fulfilled
        isConsecutiveMatch <- maxMatch >= minMatch &&
          grepl(pattern=maxSubStr, x=lowerNames2[j], fixed=TRUE) &&
          nchar(maxSubStr) > 0;

        if (isConsecutiveMatch) {
          ## merging
          mergedDf[i, "matches"] <- maxMatch
          mergedDf[i, "name2"] <- as.character(df2[j, 1]);
          #mergedDf[i, "number2"] <- df2[j, 2];

          ## don't append this row to mergedDf because already merged
          toMerge[j] <- NA;

          ## stop inner for loop here to avoid possible second match
          break;
        }
      }
    }
  }

  ## append not matched rows to mergedDf
  toMerge <- which(df2[, 1] == toMerge);

  df2 <- data.frame(name1=NA, name2=as.character(df2[toMerge, 1]), number1=NA, number2=df2[toMerge, 2], matches=0, stringsAsFactors=FALSE);

  mergedDf <- rbind(mergedDf, df2);

  return (mergedDf);
}
