## Need to work outside conda_R for this
## code from John Muschelli
# module unload conda_R
# module load R
# java_type=java-openjdk
# export _JAVA_OPTIONS="-Xms5g -Xmx6g" ## By Leo: specify some memory, otherwise java runs out of juice
# export JAVA_HOME=/usr/lib/jvm/${java_type}/jre
# export JAVA=/usr/lib/jvm/${java_type}/bin/java
# export JAVAC=/usr/lib/jvm/${java_type}/bin/javac
# export JAVAH=/usr/lib/jvm/${java_type}/bin/javah
# export JAR=/usr/lib/jvm/${java_type}/bin/jar
# export JAVA_LIBS="-L/usr/lib/jvm/${java_type}/jre/lib/amd64/server -L/usr/lib/jvm/jre-openjdk/lib/amd64/server -ljvm"
# export JAVA_CPPFLAGS="-I/usr/lib/jvm/${java_type}/jre/../include -I/usr/lib/jvm/${java_type}/jre/../include/linux"
# export JAVA_LD_LIBRARY_PATH=/usr/lib/jvm/${java_type}/jre/lib/amd64/server:/usr/lib/jvm/jre-openjdk/lib/amd64/server:/usr/lib64/lib:/usr/lib:$LD_LIBRARY_PATH
# export PATH=${PATH}:${JAVA_HOME}/../bin
# R
# install.packages('rJava')
library('xlsx')
library('sessioninfo')

load('rda/modeling_results.Rdata', verbose = TRUE)

write.xlsx(results_anova,
    file = "SupplementaryTableXX_modeling.xlsx",
    sheetName = "ANOVA model",
    append = FALSE,
    row.names = FALSE
)

write.xlsx(
    results_specificity,
    file = "SupplementaryTableXX_modeling.xlsx",
    sheetName = "Specificity model",
    append = TRUE,
    row.names = FALSE
)

write.xlsx(
    results_pairwise,
    file = "SupplementaryTableXX_modeling.xlsx",
    sheetName = "Pairwise model",
    append = TRUE,
    row.names = FALSE
)


## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────
#  setting  value
#  version  R version 3.6.1 (2019-07-05)
#  os       CentOS Linux 7 (Core)
#  system   x86_64, linux-gnu
#  ui       X11
#  language (EN)
#  collate  en_US.UTF-8
#  ctype    en_US.UTF-8
#  tz       US/Eastern
#  date     2020-02-13
#
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
#  package     * version date       lib source
#  assertthat    0.2.1   2019-03-21 [2] CRAN (R 3.6.1)
#  cli           1.1.0   2019-03-19 [2] CRAN (R 3.6.1)
#  crayon        1.3.4   2017-09-16 [2] CRAN (R 3.6.1)
#  rJava         0.9-11  2019-03-29 [1] CRAN (R 3.6.1)
#  sessioninfo * 1.1.1   2018-11-05 [1] CRAN (R 3.6.1)
#  withr         2.1.2   2018-03-15 [2] CRAN (R 3.6.1)
#  xlsx        * 0.6.2   2020-02-12 [1] CRAN (R 3.6.1)
#  xlsxjars      0.6.1   2014-08-22 [1] CRAN (R 3.6.1)
#
# [1] /users/lcollado/R/x86_64-pc-linux-gnu-library/3.6
# [2] /jhpce/shared/jhpce/core/R/3.6.1/lib64/R/library
