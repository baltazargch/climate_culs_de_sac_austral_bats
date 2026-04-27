advance_message <- function(order=NULL, species=NULL, action=NULL, appendLF = FALSE){
  stopifnot(!is.null(order), !is.null(species),
            !is.null(action), is.logical(appendLF))

  mess <- str_pad(paste0(order, ': ', species, ' (', action, ')'), width = 120, 'right')
  message('\r', mess, appendLF = appendLF)
  flush.console()
}
