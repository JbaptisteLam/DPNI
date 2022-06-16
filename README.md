# Non-Invasive Prenatal Diagnostic
Tools to estimate foetal fraction and performed non-invasive pre-natal diagnostic
# Seqff 
This repository contains also supplementaries files from seqff tool publication
Kim SK, Hannum G, Geis J, Tynan J, Hogg G, Zhao C, Jensen TJ, Mazloom AR, Oeth P, Ehrich M, van den Boom D, Deciu C. Determination of fetal DNA fraction from the plasma of pregnant women using sequence read counts. Prenat Diagn. 2015 Aug;35(8):810-5. doi: 10.1002/pd.4615. Epub 2015 Jun 3. PMID: 25967380.
> **Note:** Due to issues in args recognition arguments have been set as:
> 
|                |Arguments                    |
|----------------|-----------------------------|
| --i            |Absolute path of folder containing input files           |
| --j           modified |Name of input file (samfile without header or tabulated read                counts ordered by genomic coordinates, could be found in SupplentalTable1.csv in seqff folder)     |   
| --o             | Name of output file (it will be a csv file)
| --d             | Aboslute path of folder which will contain output
| --t             | Datatype ('sam' or 'counts' conf j arg)