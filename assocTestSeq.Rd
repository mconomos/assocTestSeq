\name{assocTestSeq}
\alias{assocTestSeq}
\title{Aggregate Association Testing with Sequencing Data}
\description{\code{assocTestSeq} performs aggregate association tests with sequencing data using the null model fit with \code{\link{fitNullMM}} or \code{\link{fitNullReg}}.}
\usage{
assocTestSeq(seqData, nullModObj, aggVarList, AF.sample = NULL,
            AF.range = c(0,1), weight.beta = c(0.5, 0.5), weight.user = NULL,
            test = "Burden", burden.test = "Score", rho = 0, 
            pval.method = "kuonen", verbose = TRUE)
}
\arguments{
    \item{seqData}{An object of class \code{SeqVarData} from the package \code{SeqVarTools} containing the sequencing genotype data for variants and samples to be used for the analysis.}
    \item{nullModObj}{A null model object returned by \code{fitNullMM} when using mixed models or \code{fitNullReg} when using linear or logistic regression.}
    \item{aggVarList}{A list specifying the variant aggregation units to be tested; each element of the list represents one aggregate test. Each element of the list should be a data.frame that contains at least two columns: variant.id matching the variant.id in seqData for the variants that should be aggregated, and allele.index specifying which alternate allele in seqData is to be tested at that variant.id.  Multiple alternate alleles can be included at the same variant location by including multiple rows with the same variant.id and different allele.index values.}
    \item{AF.sample}{A vector of sample.id values specifying which samples should be used for allele frequency calculation. When NULL (the default), all samples included in the test are used. Allele frequency calculation will affect variant inclusion based on \code{AF.range} and variant weighting based on \code{weight.beta}.}
    \item{AF.range}{A numeric vector of length two specifying the lower and upper bounds on the alternate allele frequency for variants to be included in the analysis.  Variants with alternate allele frequencies outside of this range are given a weight of 0 (i.e. excluded).}
    \item{weight.beta}{A numeric vector of length two specifying the two parameters of the Beta distribution used to determine variant weights; weights are given by \code{dbeta(AF, a, b), where AF is the alternate allele frequency, and a and b are the two parameters specified here. \code{weight.beta = c(1,25)} gives the Wu weights; \code{weight.beta} = c(0.5, 0.5)} is proportional to the Madsen-Browning weights; and \code{weight.beta = c(1,1,)} gives a weight of 1 to all variants. This input is ignored when \code{weight.user} is not NULL.}
    \item{weight.user}{A character string specifying the name of a variable in the variantData slot of the seqData object to be used as variant weights.  When left NULL (the default), the weights specified by \code{weight.beta} will be used.}
    \item{test}{A character string specifying the type of test to be performed. The possibilities are "Burden" (default) or "SKAT".  When this is set to "SKAT" and the parameter \code{rho} has multiple values, a SKAT-O test is performed.}
    \item{burden.test}{A character string specifying the type of Burden test to perform when \code{test} = "Burden".  The possibilities are "Score", "Wald", and "Firth". "Score" can be used for any \code{nullModObj}. "Wald" can not be used when the \code{nullModObj} is from a mixed model with a binary outcome variable. "Firth" can only be used when the \code{nullModObj} is from a logistic regression with a binary outcome variable.}
    \item{rho}{A numeric value (or vector of numeric values) in [0,1] specifying the rho parameter for SKAT. When rho = 0, a standard SKAT test is performed. When rho = 1, a score burden test is performed. When rho is a vector of values, SKAT-O is performed using each of those values as the search space for the optimal rho.}
    \item{pval.method}{A character string specifying which method to use to calculate SKAT p-values. "kuonen" (the default) uses a saddlepoint method; "davies" uses numerical integration; and "liu" uses a moment matching approximation.}
    \item{verbose}{Logical indicator of whether updates from the function should be printed to the console; the default is TRUE.}
}
\details{
    None yet
}
\value{A list with the following items:
    \item{param}{A list with model parameters including:}
    \item{AF.range}{The lower and upper bounds on the alternate allele frequency for variants that were included in the analysis.}
    \item{weight.beta}{The two parameters of the Beta distribution used to determine variant weights if used, NULL otherwise.}
    \item{weight.user}{A character string specifying the name of the variable in the variantData slot of the seqData object used as variant weights if used, NULL otherwise.}
    \item{family}{Either "gaussian" for a continous outcome or "binomial" for a binary outcome.}
    \item{mixedmodel}{Logical indicating whether or not a mixed model was used to fit the null model.}
    \item{test}{Specifies whether Burden, SKAT, or SKAT-O tests were performed.}
    \item{burden.test}{If test = "Burden", specifies if Score, Wald, or Firth tests were performed.}
    \item{rho}{The values of rho used in the SKAT or SKAT-O test.}
    \item{pval.method}{The p-value calculation method used in SKAT or SKAT-O tests.}

    \item{nsample}{A list with the following values:}
    \item{analysis}{The number of samples included in the analysis.}
    \item{AF}{The number of samples used to calculate allele frequencies.}

    \item{results}{A data.frame containing the results from the main analysis. Each row is a separate aggregate test:}
    \item{n.site}{The number of variant sites included in the test.}
    \item{n.sample.alt}{The number of samples with an observed alternate allele at any variant in the aggregate set.}
    If \code{test} is "Burden":
    \item{burden.skew}{The skewness of the burden value for all samples.}
    If \code{burden.test} is "Score":
    \item{Score}{The value of the score function}
    \item{Var}{The variance of the score function}
    \item{Score.stat}{The score chi-squared test statistic}
    \item{Score.pval}{The score p-value}
    If \code{burden.test} is "Wald":
    \item{Est}{The effect size estimate for a one unit increase in the burden value}
    \item{SE}{The estimated standard error of the effect size estimate}
    \item{Wald.stat}{The Wald chi-squared test statistic}
    \item{Wald.pval}{The Wald p-value}
    If \code{burden.test} is "Firth":
    \item{Est}{The effect size estimate for a one unit increase in the burden value}
    \item{SE}{The estimated standard error of the effect size estimate}
    \item{Firth.stat}{The Firth test statistic}
    \item{Firth.pval}{The Firth p-value}
    If \code{test} is "SKAT":
    \item{Q_rho}{The SKAT test statistic for the value of rho specified. There will be as many of these variables as there are rho values chosen.}
    \item{pval_rho}{The SKAT p-value for the value of rho specified.  There will be as many of these variables as there are rho values chosen.}
    \item{err_rho}{Takes value 1 if there was an error in calculating the p-value for the value of rho specified when using the "kunonen" or "davies" methods; 0 otherwise. When there is an error, the p-value returned is from the "liu" method. There will be as many of these variables as there are rho values chosen.}
    When \code{length(rho) > 1} and SKAT-O is performed:
    \item{min.pval}{The minimum p-value among the p-values calculated for each choice of rho.}
    \item{opt.rho}{The optimal rho value; i.e. the rho value that gave the minimum p-value.}
    \item{pval_SKATO}{The SKAT-O p-value after adjustment for searching across multiple rho values.}

    \item{variantInfo}{A list with as many elements as aggregate tests performed. Each element of the list is a data.frame providing information on the variants used in the aggregate test with results presented in the corresponding row of \code{results}. Each of these data.frames has the following information:}
    \item{variantID}{The variant.id value from seqData.}
    \item{allele}{The index of the allele in seqData.}
    \item{chr}{The chromosome the variant is located on.}
    \item{pos}{The position of the variant on the chromosome.}
    \item{n.obs}{The number of samples with observed genotype values at the variant.}
    \item{freq}{The allele frequency calculated using the samples specified by \code{AF.sample} (or all samples if \code{AF.sample} is NULL) of the alternate allele given by the allele index at the variant.}
    \item{weight}{The weight assigned to the variant in the analysis. A weight of 0 means the variant was excluded.}
}
%\references{}
\author{Matthew P. Conomos}
\note{None now}
%\examples{
%}
\keyword{association}
