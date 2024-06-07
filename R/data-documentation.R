#' Meta-analysis:
#'
#'
#' @name hips
#'
#' @docType data
#'
#' @description
#'
#' Meta-analysis of 15 studies investigating total hip replacement to compare the risk of revision of
#' cemented and uncemented implantfixation modalities, by pooling
#' treatment effectestimates from OS and RCTs.
#'
#' @format  A dataframe with 15 rows and 12 columns. Each row represents study results, the columns are:
#'    \describe{
#'     \item{Study}{Author and year.}
#'     \item{Study_type}{Study desing.}
#'     \item{N_of_revisions}{Number of revisions.}
#'     \item{Total_cemented}{Total number of cemmented cases.}
#'     \item{N_of_revisions_uncemented}{Number of uncemented revisions.}
#'     \item{Total_uncemented}{Total number of uncemmented cases.}
#'     \item{Relative_risks_computed}{RR calculated from the two by two table.}
#'     \item{L95CI}{Lower 95prc CI}
#'     \item{U95CI}{Upper 95prc CI}
#'     \item{mean_age}{Mean age of the study}
#'     \item{proportion_of_women}{Proportion of women in the study.}
#'     \item{Follow_up}{Time to follow-up in years.}
#'     }
#'
#' @source  Schnell-Inderst P, Iglesias CP, Arvandi M, Ciani O, Matteucci Gothe R, Peters J, Blom AW, Taylor RS and Siebert U (2017). A bias-adjusted evidence synthesis of RCT and observational data: the case of total hip replacement. Health Econ. 26(Suppl. 1): 46–69.
#'
#' @keywords datasets
#'
#'

NULL

#' Meta-analysis: Observational studies assessing the relationship of
#' a positive ICPC (Isolated Choroid Plexus Cyst) on Trisomy 21
#'
#' @name trisomy21
#'
#' @docType data
#'
#' @description
#'
#' Meta-analysis of 22 Observational Studies from PubMed, Cochrane Library and SciELO databases
#' that assessed the relationship of a positive ICPC (Isolated Choroid Plexus Cyst) on Trisomy 21
#'
#' @format  A dataframe with 22 rows and 6 columns. Each row represents study results, the columns are:
#'    \describe{
#'     \item{year}{Year of publication.}
#'     \item{author}{Principal author of the publication.}
#'     \item{y}{Number of cases of ICPC with Trisomy 21.}
#'     \item{n}{Total number o cases with ICPC.}
#'     \item{mean.GA}{Mean gestational time in weeks.}
#'     \item{study.design}{Study design: prospective or retrospective cohort.}
#'     }
#'
#' @source  Kürten C, Knippel A, Verde P, Kozlowski P. A Bayesian risk analysis for Trisomy 21 in isolated choroid plexus cyst: combining a prenatal database with a meta-analysis. J Matern Fetal Neonatal Med. 2019 Jun 11:1-9. doi: 10.1080/14767058.2019.1622666. Epub ahead of print. PMID: 31113245.
#'
#' @keywords datasets
#'
#'




NULL

#' Meta-analysis: Observational studies assessing the impact of
#' risk factors on the severity and mortality of COVID-19 cases
#'
#' @name covid19
#'
#' @docType data
#'
#' @description
#'
#' Meta-analysis of 40 Observational Studies from PubMed, Cocharane Library and SciELO databases
#' that assessed the impact of diabetes, hypertension, cardiovascular disease,
#' and the use of ACEI/ARB on severity and mortality of COVID-19 cases.
#'
#' @format  A dataframe with 89 rows and 12 columns. Each row represents study results, the columns are:
#'    \describe{
#'     \item{author}{Principal author and year of publication.}
#'     \item{endpoint}{Endoint: severity or mortality.}
#'     \item{risk.factor}{Possible risk factors: diabetes, hypertension, cardiovascular, ACE_ARB.}
#'     \item{event.e}{Number of events in the group with risk factor.}
#'     \item{n.e}{Number of patients in the group with risk factor.}
#'     \item{event.c}{Number of events in the group without risk factor.}
#'     \item{n.c}{Number of patients in the group with risk factor.}
#'     \item{design}{Study design: Case Series, Cross Sectional and Retrospective Cohort.}
#'     \item{TE}{Log Odds Ratio}
#'     \item{seTE}{Standard Error of the Log Odds Ratio}
#'     \item{logitPc}{Logit transformation of the proportion of events in the control group.}
#'     \item{N}{Total number of patients.}
#'     }
#'
#' @source  de Almeida-Pititto, B., Dualib, P.M., Zajdenverg, L. et al. Severity and mortality of COVID 19 in patients with diabetes, hypertension and cardiovascular disease: a meta-analysis. Diabetol Metab Syndr 12, 75 (2020). https://doi.org/10.1186/s13098-020-00586-4
#'
#' @keywords datasets

NULL

#' Meta-analysis of 31 randomized controled trials (RCTs) with reported discrepancies
#'
#' @name stemcells
#'
#' @docType data
#'
#' @description
#'
#' Meta-analysis of 31 randomized controled trials (RCTs) of two treatment groups
#' of heart disease patients, where the treatment group received bone marrow stem
#' cells and the control group a placebo treatment.
#'
#'
#' @format  A matrix with 31 rows and 5 columns. Each row represents study results, the columns are:
#'    \describe{
#'     \item{trial}{ID name of the trial.}
#'     \item{effect.size}{treatment effect is measured as the difference of the
#'           ejection fraction between groups, which measures the improvement of
#'           left ventricular function in the heart.}
#'     \item{se.effect}{Standard Error of the effect.size.}
#'     \item{sample.size}{Total number of patients in the trial.}
#'     \item{n.discrep}{Number of detected discrepancies in the published trial.
#'     Discrepancies are defined as two or more reported facts that cannot both be
#'     true because they are logically or mathematically incompatible.}
#'     }
#'
#' @references Verde, P. E. (2017) Two Examples of Bayesian Evidence Synthesis with the Hierarchical Meta-Regression Approach. Chap.9, pag 189-206. Bayesian Inference, ed. Tejedor, Javier Prieto. InTech.
#'
#'
#' @source  Nowbar, A N, et al. (2014) Discrepancies in autologous bone marrow stem cell trials and enhancemen of ejection fraction (DAMASCENE): weighted regression and meta-analysis. BMJ, 348,1-9.
#'
#' @keywords datasets

NULL



#' Individual participant data for diabetic patients
#'
#' @name healingipd
#'
#' @docType data
#'
#' @description
#'
#' Prospective cohort study.
#'
#' @format  A dataframe with 260 rows and 18 columns. Each row represents a patient,
#' the columns are:
#'    \describe{
#'     \item{healing.without.amp}{Outcome variable: Healing without amputation with in one year.}
#'     \item{duration_lesion_days}{Duration of leasions in days at baseline.}
#'     \item{PAD}{Peripheral arterial disease yes/no.}
#'     \item{neuropathy}{Neuropathy yes/no.}
#'     \item{first.ever.lesion}{First ever lesion yes/no.}
#'     \item{no.continuous.care}{No continuous care yes/no.}
#'     \item{male}{yes/no.}
#'     \item{diab.typ2}{Diabetes type 2 yes/no.}
#'     \item{insulin}{Insulin dependent yes/no.}
#'     \item{HOCHD}{HOCHD yes/no.}
#'     \item{HOS}{HOCHD yes/no.}
#'     \item{CRF}{CRF yes/no.}
#'     \item{dialysis}{Dialysis yes/no.}
#'     \item{DNOAP}{DNOAP yes/no.}
#'     \item{smoking.ever}{Ever smoke yes/no.}
#'     \item{age}{Age at baseline in years.}
#'     \item{diabdur}{Diabetes duration at baseline.}
#'     \item{wagner.class}{Wagner score 1-2 vs. 3-4-5.}
#'                 }
#'
#' @references Verde, P.E. (2018) The Hierarchical Meta-Regression Approach and Learning from Clinical Evidence. Technical Report.
#'
#' @source Morbach, S, et al. (2012). Long-Term Prognosis of Diabetic Foot Patients and Their Limbs: Amputation and death over the course of a decade,Diabetes Care, 35, 10, 2012-2017.
#'
#' @keywords datasets

NULL

#' Efficacy of diabetic foot healing using adjuvant treatments
#'
#' @name healing
#'
#' @docType data
#'
#' @description
#'
#'
#' Meta-analysis of 35 randomized controlled trials investigating the effectiveness
#' in the application of adjuvant therapies for diabetic patients compared to medical
#' routine care, where the endpoint was healing without amputations within a period
#' less than or equal to one year.
#'
#' @format  A matrix with 35 rows and 9 columns. Each row represents study results,
#' the columns are:
#'    \describe{
#'     \item{Study}{Name of the first author and year.}
#'     \item{n_t}{Number of patients in the treatment group.}
#'     \item{n_c}{Number of patients in the control group.}
#'     \item{y_t}{Number of heal patients in the treatment group.}
#'     \item{y_c}{Number of heal patients in the control group.}
#'     \item{ndrop}{Total number of drop out patients.}
#'     \item{fup_weeks}{Length of followup in weeks.}
#'     \item{PAD}{Inclusion of patients with peripheral arterial disease.}
#'     \item{wagner_4}{Inclusion of patients with Wagner score 3 and 4.}
#'     }
#'
#' @references Verde, P.E. (2018) The Hierarchical Meta-Regression Approach and Learning
#' from Clinical Evidence. Technical Report.
#'
#' @source  The data were obtainded from: Centre for Clinical Practice at NICE (UK and others)
#' (2011), Clinical guideline 119. Diabetic foot problems: Inpatient Management of Diabetic Foot Problems.
#'  Tech. rep., National Institute for Health and Clinical Excellence.
#'
#' @keywords datasets

NULL

#' Efficacy of Pneumococcal Polysaccharide Vaccine in Preventing Invasive Pneumococcal Disease
#'
#' @name ppvipd
#'
#' @docType data
#'
#' @description
#'
#' PPV23 (23-valent pneumococcal polysaccharide vaccine) with 3 Randomized Clinical Trials;
#' 5 Cohort Studies and 3 Case-Control Studies.
#'
#' The outcome variable IPD (Invasive Pneumococcal Disease).
#'
#' @format  A matrix with 11 rows and 6 columns. Each row represents study results, the columns are:
#'    \describe{
#'     \item{name}{Name of the first author and year.}
#'     \item{TE}{Treatment Effect as Log Odds Ratio.}
#'     \item{seTE}{Standard Error of the TE.}
#'     \item{n.v}{Number of patients in the vaccination group.}
#'     \item{n.c}{Number of patients in the control group.}
#'     \item{design}{Description of the study design.}
#'     }
#'
#' @references Falkenhorst, G., Remschmidt, C., Harder, T., Hummers-Pradier, E., Wichmann, O., and Bogdan, C. (2017) Effectiveness of the 23-Valent Pneumococcal Polysaccharide Vaccine(PPV23) against Pneumococcal Disease in the Elderly: Systematic Review and Meta-Analysis. PLoS ONE 12(1): e0169368. doi:10.1371/journal.pone.0169368.
#'
#' @references Verde, P.E. and Curcio, D. (2017) Hierarchical Meta-Regression Modelling: The Case of The Pneumococcal Polysaccharide Vaccine. Technical Report.
#'
#'
#' @source  The data were obtainded from: Falkenhorst et al. (2017).
#'
#' @keywords datasets

NULL

#' Efficacy of Pneumococcal Polysaccharide Vaccine in Preventing Community Acquired Pneumonia
#'
#' @name ppvcap
#'
#' @docType data
#'
#' @description
#'
#' PPV23 (23-valent pneumococcal polysaccharide vaccine) with 16 Randomized Clinical Trials
#' (RCTs); outcome variable CAP (community-acquired pneumonia).
#'
#' This data frame corresponds to 16 randomized control trials (RCTs) reporting efficacy of
#' the PPV (Pneumococcal Polysaccharide) vaccine in preventing CAP (community acquired pneumonia).
#' The data frame contains the evaluation of Risk of Bias (RoB) of the trials and some study population
#' characteristics.
#'
#' @format  A matrix with 16 rows and 18 columns. Each row represents study results, the columns are:
#'    \describe{
#'     \item{Name_Year}{Name of the first author and year.}
#'     \item{Year}{Year of publication.}
#'     \item{yt}{Number of infections in the intervention group.}
#'     \item{nt}{Number of patients in the intervention group.}
#'     \item{yc}{Number of infections in the control group.}
#'     \item{nc}{Number of patients in the control group.}
#'     \item{TE}{Treatment Effect as Log Odds Ratio.}
#'     \item{seTE}{Standard Error of the TE.}
#'     \item{logitPc}{Observed baseline rate in logit scale.}
#'     \item{N}{Total sample size.}
#'     \item{Study_Design}{Description of the study design.}
#'     \item{Intervention}{Type of vaccine used for itervention.}
#'     \item{Valency}{0 = PPV23; 1 = PPV-Other.}
#'     \item{low_income}{Indicates low income patients population with 0 = no; 1 = yes.}
#'     \item{R1}{Random sequence generation (selection bias: low;high;unclear.}
#'     \item{R2}{Allocation concealment (selection bias): low;high;unclear.}
#'     \item{R3}{Confounding: low;high;unclear.}
#'     \item{R4}{Blinding of participants and personnel (performace bias): low;high;unclear.}
#'     \item{R5}{Blinding of outcome assessment (detection bias): low;high;unclear.}
#'     \item{R6}{Incomplete outcome data (attrition bias): low;high;unclear.}
#'     \item{R7}{Selective reporting (reporting bias): low;high;unclear.}
#'     \item{Participants}{Comments on patients characteristics.}
#'     }
#'
#' @references Moberley, S., Holden, J., Tatham, D., and Andrews, R. (2013), Vaccines for preventing pneumococcal infection in adults., Cochrane Database of Systematic Reviews, Issue 1. Art. No.: CD000422. DOI:10.1002/14651858.CD000422.pub3.
#'
#' @references Verde, P.E. and Curcio, D. (2017) Hierarchical Meta-Regression Modelling: The Case of The Pneumococcal Polysaccharide Vaccine. Technical Report.
#'
#' @source  The data were obtainded from: Moberley et al. (2013).
#'
#' @keywords datasets
NULL

