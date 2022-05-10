curve.one <- function(x, C, SD){
    return(
        (C*exp(-(x-0)**2/(2*SD**2)))
    )
}

curve.two <- function(x, C, SD){
    return(
        (C*exp(-(x-0)**2/(2*SD**2)))
    )
}

curve.three <- function(x, C, SD){
    return(
        (C*exp(-(x-0)**2/(2*SD**2)))
    )
}

all.curves <- function(x, nr.of.curves, C1, C2, C3=NULL, SD1, SD2, SD3=NULL){
    if(nr.of.curves == 2){
        return(
            (C1*exp(-(x-0)**2/(2*SD1**2))+
             C2*exp(-(x-0)**2/(2*SD2**2))
            )
        )
    } else if(nr.of.curves == 3){
        return(
            (C1*exp(-(x-0)**2/(2*SD1**2))+
             C2*exp(-(x-0)**2/(2*SD2**2))+
             C3*exp(-(x-0)**2/(2*SD3**2))
            )
        )
    }
}