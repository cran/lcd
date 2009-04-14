# Some utility function in the algorithm

setClass("sep.pair",
         representation(u="character", v="character", s="character"),
         package = "lcd")

setMethod("show", signature("sep.pair"),
          function(object) print(c(u=object@u, v=object@v, s=object@s))
          )

setMethod("all.equal", signature(target="sep.pair", current="sep.pair"),
          function(target, current) {
              if(target@u == current@u) {
                  return (target@v == current@v)
              } else {
                  return (target@u == current@v && target@v == current@u)
              }
          })




