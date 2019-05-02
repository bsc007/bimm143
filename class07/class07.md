class07
================
Benjamin Cho
April 23, 2019

fUNCTIONS REVISTED
==================

``` r
source("http://tinyurl.com/rescale-R")
```

``` r
rescale(1:10)
```

    ##  [1] 0.0000000 0.1111111 0.2222222 0.3333333 0.4444444 0.5555556 0.6666667
    ##  [8] 0.7777778 0.8888889 1.0000000

Try the rescale2() function that catches string inputs

``` r
rescale2(c(1:10), "string")
```

    ##  [1] 0.0000000 0.1111111 0.2222222 0.3333333 0.4444444 0.5555556 0.6666667
    ##  [8] 0.7777778 0.8888889 1.0000000

``` r
x <- c( 1, 2, NA, 3, NA)
y <- c(NA, 3, NA, 3, 4)
```

``` r
is.na(x)
```

    ## [1] FALSE FALSE  TRUE FALSE  TRUE

``` r
is.na(y)
```

    ## [1]  TRUE FALSE  TRUE FALSE FALSE

``` r
is.na(x) & is.na(y)
```

    ## [1] FALSE FALSE  TRUE FALSE FALSE

Take the sum() to find out how many TRUE values we have and thus how many NAs we had in both x and y

``` r
sum(is.na(x) & is.na(y))
```

    ## [1] 1

Now I can make this into our first function

``` r
both_na <- function(x,y) {
  sum(is.na(x) & is.na(y))
}
```

``` r
both_na(x,c(NA,3,NA,2, NA))
```

    ## [1] 2

Test

``` r
x <- c(NA, NA, NA)
y1 <- c( 1, NA, NA)
y2 <- c( 1, NA, NA, NA)
```

``` r
both_na(x,y2)
```

    ## Warning in is.na(x) & is.na(y): longer object length is not a multiple of
    ## shorter object length

    ## [1] 3

``` r
#it matches the shorter object to make it the same length of the longer one. It will recycle back the shorter object
```

``` r
y3 <- c(1, NA, NA, NA, NA, NA)
both_na(x,y3)
```

    ## [1] 5

``` r
#tell me if they are equal
3==3
```

    ## [1] TRUE

``` r
3>2
```

    ## [1] TRUE

``` r
# not equal
3 != 2
```

    ## [1] TRUE

``` r
length(x)
```

    ## [1] 3

``` r
length(x) != length(y2)
```

    ## [1] TRUE

now lets try the both\_na2 function because it was a different length

``` r
#both_na2(x,y2)
```

``` r
which(c(F,F,T,F, T))
```

    ## [1] 3 5

``` r
#which parts of the vector are true
which(is.na(c(1,2,NA, 4)))
```

    ## [1] 3

``` r
x <- c( 1, 2, NA, 3, NA)
y <- c(NA, 3, NA, 3, NA)

both_na3(x,y)
```

    ## Found 2 NA's at position(s):3, 5

    ## $number
    ## [1] 2
    ## 
    ## $which
    ## [1] 3 5

Intersect functions
===================

``` r
df1
```

    ##     IDs exp
    ## 1 gene1   2
    ## 2 gene2   1
    ## 3 gene3   1

``` r
df2
```

    ##     IDs exp
    ## 1 gene2  -2
    ## 2 gene4  NA
    ## 3 gene3   1
    ## 4 gene5   2

``` r
x <- df1$IDs
y <- df2$IDs

x
```

    ## [1] "gene1" "gene2" "gene3"

``` r
y
```

    ## [1] "gene2" "gene4" "gene3" "gene5"

``` r
?intersect
```

    ## starting httpd help server ... done

``` r
intersect(x,y)
```

    ## [1] "gene2" "gene3"

``` r
which(x %in% y)
```

    ## [1] 2 3

``` r
y[y%in%x]
```

    ## [1] "gene2" "gene3"

A working snippet

Use Rstudio shortcut CODE &gt; EXTRACT FUNCTION to turn our code into a function

``` r
gene_intersect <- function(x, y) {
  cbind( x[ x %in% y ], y[ y %in% x ] )
}
```

``` r
gene_intersect(df1$exp,df1$exp)
```

    ##      [,1] [,2]
    ## [1,]    2    2
    ## [2,]    1    1
    ## [3,]    1    1

gene\_intersect2(df1,df2)

``` r
gene_intersect2(df1,df2)
```

    ##     IDs exp df2[df2$IDs %in% df1$IDs, "exp"]
    ## 2 gene2   1                               -2
    ## 3 gene3   1                                1

``` r
gene_intersect3(df1,df2)
```

    ##     IDs exp exp2
    ## 2 gene2   1   -2
    ## 3 gene3   1    1

``` r
merge(df1,df2,by="IDs")
```

    ##     IDs exp.x exp.y
    ## 1 gene2     1    -2
    ## 2 gene3     1     1

Find the average grade dropping the worst homework score

``` r
grade <- function(x) {
  exclude<-x[x != min(x)]
  mean(exclude)
}
```

``` r
bob<- c(100,100,100,100,100, 100, 100, 90)
grade(bob)
```

    ## [1] 100

You install.package a package once, but you have to use library to be able to retrieve the function
