#' Title
#'
#' @param icog_onco_score_infor_one
#' @param second.num
#'
#' @return
#' @export
#'
#' @examples
MetaPfunction <- function(icog_onco_score_infor_one,second.num){
  score.icog <- as.numeric(icog_onco_score_infor_one[1:(second.num)])
  infor.icog <- matrix(as.numeric(icog_onco_score_infor_one[(second.num+1):
                                                              (second.num+second.num^2) ]),
                       ncol = second.num)
  start <- second.num+second.num^2
  score.onco <- as.numeric(icog_onco_score_infor_one[(1+start):
                                                       (second.num+start)])
  infor.onco <- matrix(as.numeric(icog_onco_score_infor_one[(second.num+1+start):
                                                              (second.num+second.num^2+start) ]),ncol=second.num)


  meta.result <- ScoreMetaAnalysis(score.icog,infor.icog,
                                   score.onco,infor.onco)
  score.meta <- t(meta.result[[1]])
  infor.meta <- meta.result[[2]]
  DisplayFixedScoreTestResult(score.meta,infor.meta)
}
