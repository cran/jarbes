#' Meta-analysis: Long-COVID Health Outcomes
#'
#' @name longcovid
#'
#' @docType data
#'
#' @description
#' This dataset is based on a comprehensive meta-analysis of 33 studies, sourced from various databases,
#' including the Cochrane COVID-19 Study Register (comprising the Cochrane Central Register of Controlled Trials,
#' Medline, Embase, clinicaltrials.gov, the World Health Organization's International Clinical Trials Registry Platform,
#' and medRxiv) and the World Health Organization’s COVID-19 research database. The analysis focused on evaluating
#' health outcomes related to Long-COVID in controlled studies. Specifically, it examines the health outcomes in
#' terms of incident medicinal diagnoses.
#'
#' The dataset includes the assessment of risk of bias based on the Joanna Briggs Institute (JBI) tool for cohort studies,
#' along with various participant and study details such as sample size, effect type, follow-up time, and disease severity.
#'
#' @format A data frame with 271 rows and 27 columns. Each row represents the results of a single study.
#' The columns include:
#' \describe{
#'   \item{study}{Name of the first author and publication year.}
#'   \item{category}{Category of the health outcome.}
#'   \item{outcome_disease}{Definition of the health outcome or disease.}
#'   \item{data_source}{Type of data source: Administrative data, Health records, Patients claims, Survey, Combination of health records and claims. }
#'   \item{sample_size}{Total number of participants.}
#'   \item{effect_type}{Type of effect reported: RR (Relative Risk), HR (Hazard Ratio), or OR (Odds Ratio).}
#'   \item{effect}{Estimated effect based on the effect type.}
#'   \item{TE}{Logarithm of the estimated effect.}
#'   \item{seTE}{Standard error of the logarithm of the estimated effect.}
#'   \item{rate_control}{Event rate in the control group.}
#'   \item{follow_up_time}{Follow-up time in weeks.}
#'   \item{mean_age}{Mean age of the participants.}
#'   \item{disease_severity}{Indicator for inclusion of severe or critical disease participants ("no" or "yes").}
#'   \item{reinfection}{Indicator for inclusion of reinfected participants ("no" or "yes").}
#'   \item{no_of_confounders}{Number of confounders for which adjustments were made in the study.}
#'   \item{uncertainty_of_confounders}{high if ROB4 OR ROB5 is high or unclear or low otherwise.}
#'   \item{list_of_confounders}{List of confounders considered in the study.}
#'   \item{ROB1}{Were the two groups similar and recruited from the same population?}
#'    \item{ROB2}{Were the exposures measured similarly to assign participants to exposed and unexposed groups?}
#'    \item{ROB3}{Was the exposure measured in a valid and reliable way?}
#'    \item{ROB4}{Were confounding factors identified?}
#'    \item{ROB5}{Were strategies to address confounding factors stated?}
#'    \item{ROB6}{Were the groups/participants free of the outcome at the start of the study (or at the moment of exposure)?}
#'    \item{ROB7}{Were the outcomes measured in a valid and reliable way?}
#'    \item{ROB8}{Was the follow-up time reported and sufficient to allow outcomes to occur?}
#'    \item{ROB9}{Was follow-up complete, and if not, were the reasons for loss to follow-up described and explored?}
#'    \item{ROB10}{Were strategies to address incomplete follow-up utilized?}
#'    \item{ROB11}{Was appropriate statistical analysis used?}
#' }
#'
#' @source Franco JVA, Garegnani LI, Metzendorf MI, Heldt K, Mumm R, Scheidt-Nave C.
#' Post-COVID-19 conditions in adults: systematic review and meta-analysis of health outcomes
#' in controlled studies. BMJ Medicine. 2024;3:e000723.
#'

#' @keywords datasets
#'
NULL


#' Individual Participant Data: Diabetic Eyes Data
#'
#'
#' @name diabetes_eyes
#'
#' @docType data
#'
#' @description
#'
#' A dataset containing detailed measurements from a study investigating the relationship between diabetes and eye health. The dataset includes patient demographics, visual acuity, and extensive macular metrics derived from optical coherence tomography (OCT) imaging.
#'
#' @format  A dataframe with 270 columns and 97 rows, where each row represents a patient. The columns include:
#'    \describe{
#'     \item{pat}{Patient ID.}
#'     \item{diabetes_type}{Indicator for diabetes (2 = diabetic (type 2), 0 = healthy).}
#'     \item{sex}{Gender of the patient (1 = Male, 2 = Female).}
#'     \item{age}{Age of the patient (years).}
#'     \item{smoker}{Smoking status (1 = Smoker, 0 = Non-smoker).}
#'     \item{weight}{Weight of the patient (kg).}
#'     \item{height}{Height of the patient (m).}
#'     \item{BMI}{Body Mass Index}
#'     \item{VISUAL_ACUITY_RIGHT_EYE}{Visual acuity for the right eye.}
#'     \item{VISUAL_ACUITY_LEFT_EYE}{Visual acuity for the left eye.}
#'     \item{CONTRAST_SENSITIVITY_RIGHT_EYE}{Measure of contrast sensitivity for the left eye.}
#'     \item{CONTRAST_SENSITIVITY_LEFT_EYE}{Measure of contrast sensitivity for the left eye.}
#'     \item{R_PAP_RNFL_N}{Measurement of the right eye, Papilla, Retinal Nerve Fiber Layer, Nasal Parafovea}
#'     \item{R_PAP_RNFL_NI}{Measurement of the right eye, Papilla, Retinal Nerve Fiber Layer, Nasal Inferior Parafovea}
#'     \item{R_PAP_RNFL_TI}{Measurement of the right eye, Papilla, Retinal Nerve Fiber Layer, Temporal Inferior Parafovea}
#'     \item{R_PAP_RNFL_T}{Measurement of the right eye, Papilla, Retinal Nerve Fiber Layer, Temporal Parafovea}
#'     \item{R_PAP_RNFL_TS}{Measurement of the right eye, Papilla, Retinal Nerve Fiber Layer, Temporal Superior Parafovea}
#'     \item{R_PAP_RNFL_G}{Measurement of the right eye, Papilla, Retinal Nerve Fiber Layer, Global Layer}
#'     \item{L_PAP_RNFL_NS}{Measurement of the left eye, Papilla, Retinal Nerve Fiber Layer, Nasal Superior Parafovea}
#'     \item{L_PAP_RNFL_N}{Measurement of the left eye, Papilla, Retinal Nerve Fiber Layer, Nasal Parafovea}
#'     \item{L_PAP_RNFL_NI}{Measurement of the left eye, Papilla, Retinal Nerve Fiber Layer, Nasal Inferior Parafovea}
#'     \item{L_PAP_RNFL_TI}{Measurement of the left eye, Papilla, Retinal Nerve Fiber Layer, Temporal Inferior Parafovea}
#'     \item{L_PAP_RNFL_T}{Measurement of the left eye, Papilla, Retinal Nerve Fiber Layer, Temporal Parafovea}
#'     \item{L_PAP_RNFL_TS}{Measurement of the left eye, Papilla, Retinal Nerve Fiber Layer, Temporal Superior Parafovea}
#'     \item{L_PAP_RNFL_G}{Measurement of the left eye, Papilla, Retinal Nerve Fiber Layer, Global Layer}
#'     \item{R_PAP_FULL_NS}{Measurement of the right eye, Papilla, Complete Retinal Thickness, Nasal Superior Parafovea}
#'     \item{R_PAP_FULL_N}{Measurement of the right eye, Papilla, Complete Retinal Thickness, Nasal Parafovea}
#'     \item{R_PAP_FULL_NI}{Measurement of the right eye, Papilla, Complete Retinal Thickness, Nasal Inferior Parafovea}
#'     \item{R_PAP_FULL_TI}{Measurement of the right eye, Papilla, Complete Retinal Thickness, Temporal Inferior Parafovea}
#'     \item{R_PAP_FULL_T}{Measurement of the right eye, Papilla, Complete Retinal Thickness, Temporal Parafovea}
#'     \item{R_PAP_FULL_TS}{Measurement of the right eye, Papilla, Complete Retinal Thickness, Temporal Superior Parafovea}
#'     \item{R_PAP_FULL_G}{Measurement of the right eye, Papilla, Complete Retinal Thickness, Global Layer}
#'     \item{L_PAP_FULL_NS}{Measurement of the left eye, Papilla, Complete Retinal Thickness, Nasal Superior Parafovea}
#'     \item{L_PAP_FULL_N}{Measurement of the left eye, Papilla, Complete Retinal Thickness, Nasal Parafovea}
#'     \item{L_PAP_FULL_NI}{Measurement of the left eye, Papilla, Complete Retinal Thickness, Nasal Inferior Parafovea}
#'     \item{L_PAP_FULL_TI}{Measurement of the left eye, Papilla, Complete Retinal Thickness, Temporal Inferior Parafovea}
#'     \item{L_PAP_FULL_T}{Measurement of the left eye, Papilla, Complete Retinal Thickness, Temporal Parafovea}
#'     \item{L_PAP_FULL_TS}{Measurement of the left eye, Papilla, Complete Retinal Thickness, Temporal Superior Parafovea}
#'     \item{L_PAP_FULL_G}{Measurement of the left eye, Papilla, Complete Retinal Thickness, Global Layer}
#'     \item{R_PAP_GCLIPL_NS}{Measurement of the right eye, Papilla, Ganglion Cell Layer and Inner Plexiform Layer combined, Nasal Superior Parafovea}
#'     \item{R_PAP_GCLIPL_N}{Measurement of the right eye, Papilla, Ganglion Cell Layer and Inner Plexiform Layer combined, Nasal Parafovea}
#'     \item{R_PAP_GCLIPL_NI}{Measurement of the right eye, Papilla, Ganglion Cell Layer and Inner Plexiform Layer combined, Nasal Inferior Parafovea}
#'     \item{R_PAP_GCLIPL_TI}{Measurement of the right eye, Papilla, Ganglion Cell Layer and Inner Plexiform Layer combined, Temporal Inferior Parafovea}
#'     \item{R_PAP_GCLIPL_T}{Measurement of the right eye, Papilla, Ganglion Cell Layer and Inner Plexiform Layer combined, Temporal Parafovea}
#'     \item{R_PAP_GCLIPL_TS}{Measurement of the right eye, Papilla, Ganglion Cell Layer and Inner Plexiform Layer combined, Temporal Superior Parafovea}
#'     \item{R_PAP_GCLIPL_G}{Measurement of the right eye, Papilla, Ganglion Cell Layer and Inner Plexiform Layer combined, Global Layer}
#'     \item{L_PAP_GCLIPL_NS}{Measurement of the left eye, Papilla, Ganglion Cell Layer and Inner Plexiform Layer combined, Nasal Superior Parafovea}
#'     \item{L_PAP_GCLIPL_N}{Measurement of the left eye, Papilla, Ganglion Cell Layer and Inner Plexiform Layer combined, Nasal Parafovea}
#'     \item{L_PAP_GCLIPL_NI}{Measurement of the left eye, Papilla, Ganglion Cell Layer and Inner Plexiform Layer combined, Nasal Inferior Parafovea}
#'     \item{L_PAP_GCLIPL_TI}{Measurement of the left eye, Papilla, Ganglion Cell Layer and Inner Plexiform Layer combined, Temporal Inferior Parafovea}
#'     \item{L_PAP_GCLIPL_T}{Measurement of the left eye, Papilla, Ganglion Cell Layer and Inner Plexiform Layer combined, Temporal Parafovea}
#'     \item{L_PAP_GCLIPL_TS}{Measurement of the left eye, Papilla, Ganglion Cell Layer and Inner Plexiform Layer combined, Temporal Superior Parafovea}
#'     \item{L_PAP_GCLIPL_G}{Measurement of the left eye, Papilla, Ganglion Cell Layer and Inner Plexiform Layer combined, Global Layer}
#'     \item{R_PAP_INLOPL_NS}{Measurement of the right eye, Papilla, Inner Nuclear Layer and Outer Plexiform Layer combined, Nasal Superior Parafovea}
#'     \item{R_PAP_INLOPL_N}{Measurement of the right eye, Papilla, Inner Nuclear Layer and Outer Plexiform Layer combined, Nasal Parafovea}
#'     \item{R_PAP_INLOPL_NI}{Measurement of the right eye, Papilla, Inner Nuclear Layer and Outer Plexiform Layer combined, Nasal Inferior Parafovea}
#'     \item{R_PAP_INLOPL_TI}{Measurement of the right eye, Papilla, Inner Nuclear Layer and Outer Plexiform Layer combined, Temporal Inferior Parafovea}
#'     \item{R_INLOPL_T}{Measurement of the right eye, Inner Nuclear Layer and Outer Plexiform Layer combined, Temporal Parafovea}
#'     \item{R_PAP_INLOPL_TS}{Measurement of the right eye, Papilla, Inner Nuclear Layer and Outer Plexiform Layer combined, Temporal Superior Parafovea}
#'     \item{R_PAP_INLOPL_G}{Measurement of the right eye, Papilla, Inner Nuclear Layer and Outer Plexiform Layer combined, Global Layer}
#'     \item{L_PAP_INLOPL_NS}{Measurement of the left eye, Papilla, Inner Nuclear Layer and Outer Plexiform Layer combined, Nasal Superior Parafovea}
#'     \item{L_PAP_INLOPL_N}{Measurement of the left eye, Papilla, Inner Nuclear Layer and Outer Plexiform Layer combined, Nasal Parafovea}
#'     \item{L_PAP_INLOPL_NI}{Measurement of the left eye, Papilla, Inner Nuclear Layer and Outer Plexiform Layer combined, Nasal Inferior Parafovea}
#'     \item{L_PAP_INLOPL_TI}{Measurement of the left eye, Papilla, Inner Nuclear Layer and Outer Plexiform Layer combined, Temporal Inferior Parafovea}
#'     \item{L_PAP_INLOPL_T}{Measurement of the left eye, Papilla, Inner Nuclear Layer and Outer Plexiform Layer combined, Temporal Parafovea}
#'     \item{L_PAP_INLOPL_TS}{Measurement of the left eye, Papilla, Inner Nuclear Layer and Outer Plexiform Layer combined, Temporal Superior Parafovea}
#'     \item{L_PAP_INLOPL_G}{Measurement of the left eye, Papilla, Inner Nuclear Layer and Outer Plexiform Layer combined, Global Layer}
#'     \item{R_PAP_ONLFIS_NS}{Measurement of the right eye, Papilla, Outer Nuclear Layer and Photoreceptor Inner Segment Layer combined, Nasal Superior Parafovea}
#'     \item{R_PAP_ONLFIS_N}{Measurement of the right eye, Papilla, Outer Nuclear Layer and Photoreceptor Inner Segment Layer combined, Nasal Parafovea}
#'     \item{R_PAP_ONLFIS_NI}{Measurement of the right eye, Papilla, Outer Nuclear Layer and Photoreceptor Inner Segment Layer combined, Nasal Inferior Parafovea}
#'     \item{R_PAP_ONLFIS_TI}{Measurement of the right eye, Papilla, Outer Nuclear Layer and Photoreceptor Inner Segment Layer combined, Temporal Inferior Parafovea}
#'     \item{R_PAP_ONLFIS_T}{Measurement of the right eye, Papilla, Outer Nuclear Layer and Photoreceptor Inner Segment Layer combined, Temporal Parafovea}
#'     \item{R_PAP_ONLFIS_TS}{Measurement of the right eye, Papilla, Outer Nuclear Layer and Photoreceptor Inner Segment Layer combined, Temporal Superior Parafovea}
#'     \item{R_PAP_ONLFIS_G}{Measurement of the right eye, Papilla, Outer Nuclear Layer and Photoreceptor Inner Segment Layer combined, Global Layer}
#'     \item{L_PAP_ONLFIS_NS}{Measurement of the left eye, Papilla, Outer Nuclear Layer and Photoreceptor Inner Segment Layer combined, Nasal Superior Parafovea}
#'     \item{L_PAP_ONLFIS_N}{Measurement of the left eye, Papilla, Outer Nuclear Layer and Photoreceptor Inner Segment Layer combined, Nasal Parafovea}
#'     \item{L_PAP_ONLFIS_NI}{Measurement of the left eye, Papilla, Outer Nuclear Layer and Photoreceptor Inner Segment Layer combined, Nasal Inferior Parafovea}
#'     \item{L_PAP_ONLFIS_TI}{Measurement of the left eye, Papilla, Outer Nuclear Layer and Photoreceptor Inner Segment Layer combined, Temporal Inferior Parafovea}
#'     \item{L_PAP_ONLFIS_T}{Measurement of the left eye, Papilla, Outer Nuclear Layer and Photoreceptor Inner Segment Layer combined, Temporal Parafovea}
#'     \item{L_PAP_ONLFIS_TS}{Measurement of the left eye, Papilla, Outer Nuclear Layer and Photoreceptor Inner Segment Layer combined, Temporal Superior Parafovea}
#'     \item{L_PAP_ONLFIS_G}{Measurement of the left eye, Papilla, Outer Nuclear Layer and Photoreceptor Inner Segment Layer combined, Global Layer}
#'     \item{R_PAP_FBBM_NS}{Measurement of the right eye, Papilla, Photoreceptor Fiber Layer and Basal Membrane combined, Nasal Superior Parafovea}
#'     \item{R_PAP_FBBM_N}{Measurement of the right eye, Papilla, Photoreceptor Fiber Layer and Basal Membrane combined, Nasal Parafovea}
#'     \item{R_PAP_FBBM_NI}{Measurement of the right eye, Papilla, Photoreceptor Fiber Layer and Basal Membrane combined, Nasal Inferior Parafovea}
#'     \item{R_PAP_FBBM_TI}{Measurement of the right eye, Papilla, Photoreceptor Fiber Layer and Basal Membrane combined, Temporal Inferior Parafovea}
#'     \item{R_PAP_FBBM_T}{Measurement of the right eye, Papilla, Photoreceptor Fiber Layer and Basal Membrane combined, Temporal Parafovea}
#'     \item{R_PAP_FBBM_TS}{Measurement of the right eye, Papilla, Photoreceptor Fiber Layer and Basal Membrane combined, Temporal Superior Parafovea}
#'     \item{R_PAP_FBBM_G}{Measurement of the right eye, Papilla, Photoreceptor Fiber Layer and Basal Membrane combined, Global Layer}
#'     \item{L_PAP_FBBM_NS}{Measurement of the left eye, Papilla, Photoreceptor Fiber Layer and Basal Membrane combined, Nasal Superior Parafovea}
#'     \item{L_PAP_FBBM_N}{Measurement of the left eye, Papilla, Photoreceptor Fiber Layer and Basal Membrane combined, Nasal Parafovea}
#'     \item{L_PAP_FBBM_NI}{Measurement of the left eye, Papilla, Photoreceptor Fiber Layer and Basal Membrane combined, Nasal Inferior Parafovea}
#'     \item{L_PAP_FBBM_TI}{Measurement of the left eye, Papilla, Photoreceptor Fiber Layer and Basal Membrane combined, Temporal Inferior Parafovea}
#'     \item{L_PAP_FBBM_T}{Measurement of the left eye, Papilla, Photoreceptor Fiber Layer and Basal Membrane combined, Temporal Parafovea}
#'     \item{L_PAP_FBBM_TS}{Measurement of the left eye, Papilla, Photoreceptor Fiber Layer and Basal Membrane combined, Temporal Superior Parafovea}
#'     \item{L_PAP_FBBM_G}{Measurement of the left eye, Papilla, Photoreceptor Fiber Layer and Basal Membrane combined, Global Layer}
#'     \item{M_R_MAC_FULL_N2}{Manual measurement of the right eye, Macula, Complete Retinal Thickness, Nasal Outer Parafovea}
#'     \item{M_R_MAC_GCLIPL_N2}{Manual measurement of the right eye, Macula, Ganglion Cell Layer and Inner Plexiform Layer combined, Nasal Outer Parafovea}
#'     \item{M_R_MAC_INLOPL_N2}{Manual measurement of the right eye, Macula, Inner Nuclear Layer and Outer Plexiform Layer combined, Nasal Outer Parafovea}
#'     \item{M_R_MAC_ONLFIS_N2}{Manual measurement of the right eye, Macula, Outer Nuclear Layer and Photoreceptor Inner Segment Layer combined, Nasal Outer Parafovea}
#'     \item{M_R_MAC_FBBM_N2}{Manual measurement of the right eye, Macula, Photoreceptor Fiber Layer and Basal Membrane combined, Nasal Outer Parafovea}
#'     \item{M_R_MAC_RNFL_N2}{Manual measurement of the right eye, Macula, Retinal Nerve Fiber Layer, Nasal Outer Parafovea}
#'     \item{M_L_MAC_FULL_N2}{Manual measurement of the left eye, Macula, Complete Retinal Thickness, Nasal Outer Parafovea}
#'     \item{M_L_MAC_GCLIPL_N2}{Manual measurement of the left eye, Macula, Ganglion Cell Layer and Inner Plexiform Layer combined, Nasal Outer Parafovea}
#'     \item{M_L_MAC_INLOPL_N2}{Manual measurement of the left eye, Macula, Inner Nuclear Layer and Outer Plexiform Layer combined, Nasal Outer Parafovea}
#'     \item{M_L_MAC_ONLFIS_N2}{Manual measurement of the left eye, Macula, Outer Nuclear Layer and Photoreceptor Inner Segment Layer combined, Nasal Outer Parafovea}
#'     \item{M_L_MAC_FBBM_N2}{Manual measurement of the left eye, Macula, Photoreceptor Fiber Layer and Basal Membrane combined, Nasal Outer Parafovea}
#'     \item{M_L_MAC_RNFL_N2}{Manual measurement of the left eye, Macula, Retinal Nerve Fiber Layer, Nasal Outer Parafovea}
#'     \item{R_MAC_FULL_S1}{Measurement of the right eye, Macula, Complete Retinal Thickness, Superior Inner Parafovea}
#'     \item{R_MAC_FULL_S2}{Measurement of the right eye, Macula, Complete Retinal Thickness, Superior Outer Parafovea}
#'     \item{R_MAC_FULL_N1}{Measurement of the right eye, Macula, Complete Retinal Thickness, Nasal Inner Parafovea}
#'     \item{R_MAC_FULL_N2}{Measurement of the right eye, Macula, Complete Retinal Thickness, Nasal Outer Parafovea}
#'     \item{R_MAC_FULL_I1}{Measurement of the right eye, Macula, Complete Retinal Thickness, Inferior Inner Parafovea}
#'     \item{R_MAC_FULL_I2}{Measurement of the right eye, Macula, Complete Retinal Thickness, Inferior Outer Parafovea}
#'     \item{R_MAC_FULL_T1}{Measurement of the right eye, Macula, Complete Retinal Thickness, Temporal Inner Parafovea}
#'     \item{R_MAC_FULL_T2}{Measurement of the right eye, Macula, Complete Retinal Thickness, Temporal Outer Parafovea}
#'     \item{R_MAC_FULL_C}{Measurement of the right eye, Macula, Complete Retinal Thickness, Center Fovea}
#'     \item{R_MAC_RNFL_S1}{Measurement of the right eye, Macula, Retinal Nerve Fiber Layer, Superior Inner Parafovea}
#'     \item{R_MAC_RNFL_S2}{Measurement of the right eye, Macula, Retinal Nerve Fiber Layer, Superior Outer Parafovea}
#'     \item{R_MAC_RNFL_N1}{Measurement of the right eye, Macula, Retinal Nerve Fiber Layer, Nasal Inner Parafovea}
#'     \item{R_MAC_RNFL_N2}{Measurement of the right eye, Macula, Retinal Nerve Fiber Layer, Nasal Outer Parafovea}
#'     \item{R_MAC_RNFL_I1}{Measurement of the right eye, Macula, Retinal Nerve Fiber Layer, Inferior Inner Parafovea}
#'     \item{R_MAC_RNFL_I2}{Measurement of the right eye, Macula, Retinal Nerve Fiber Layer, Inferior Outer Parafovea}
#'     \item{R_MAC_RNFL_T1}{Measurement of the right eye, Macula, Retinal Nerve Fiber Layer, Temporal Inner Parafovea}
#'     \item{R_MAC_RNFL_T2}{Measurement of the right eye, Macula, Retinal Nerve Fiber Layer, Temporal Outer Parafovea}
#'     \item{R_MAC_RNFL_C}{Measurement of the right eye, Macula, Retinal Nerve Fiber Layer, Center Fovea}
#'     \item{R_MAC_GCL_S1}{Measurement of the right eye, Macula, Ganglion Cell Layer, Superior Inner Parafovea}
#'     \item{R_MAC_GCL_S2}{Measurement of the right eye, Macula, Ganglion Cell Layer, Superior Outer Parafovea}
#'     \item{R_MAC_GCL_N1}{Measurement of the right eye, Macula, Ganglion Cell Layer, Nasal Inner Parafovea}
#'     \item{R_MAC_GCL_N2}{Measurement of the right eye, Macula, Ganglion Cell Layer, Nasal Outer Parafovea}
#'     \item{R_MAC_GCL_I1}{Measurement of the right eye, Macula, Ganglion Cell Layer, Inferior Inner Parafovea}
#'     \item{R_MAC_GCL_I2}{Measurement of the right eye, Macula, Ganglion Cell Layer, Inferior Outer Parafovea}
#'     \item{R_MAC_GCL_T1}{Measurement of the right eye, Macula, Ganglion Cell Layer, Temporal Inner Parafovea}
#'     \item{R_MAC_GCL_T2}{Measurement of the right eye, Macula, Ganglion Cell Layer, Temporal Outer Parafovea}
#'     \item{R_MAC_GCL_C}{Measurement of the right eye, Macula, Ganglion Cell Layer, Center Fovea}
#'     \item{R_MAC_IPL_S1}{Measurement of the right eye, Macula, Inner Plexiform Layer, Superior Inner Parafovea}
#'     \item{R_MAC_IPL_S2}{Measurement of the right eye, Macula, Inner Plexiform Layer, Superior Outer Parafovea}
#'     \item{R_MAC_IPL_N1}{Measurement of the right eye, Macula, Inner Plexiform Layer, Nasal Inner Parafovea}
#'     \item{R_MAC_IPL_N2}{Measurement of the right eye, Macula, Inner Plexiform Layer, Nasal Outer Parafovea}
#'     \item{R_MAC_IPL_I1}{Measurement of the right eye, Macula, Inner Plexiform Layer, Inferior Inner Parafovea}
#'     \item{R_MAC_IPL_I2}{Measurement of the right eye, Macula, Inner Plexiform Layer, Inferior Outer Parafovea}
#'     \item{R_MAC_IPL_T1}{Measurement of the right eye, Macula, Inner Plexiform Layer, Temporal Inner Parafovea}
#'     \item{R_MAC_IPL_T2}{Measurement of the right eye, Macula, Inner Plexiform Layer, Temporal Outer Parafovea}
#'     \item{R_MAC_IPL_C}{Measurement of the right eye, Macula, Inner Plexiform Layer, Center Fovea}
#'     \item{R_MAC_INL_S1}{Measurement of the right eye, Macula, Inner Nuclear Layer, Superior Inner Parafovea}
#'     \item{R_MAC_INL_S2}{Measurement of the right eye, Macula, Inner Nuclear Layer, Superior Outer Parafovea}
#'     \item{R_MAC_INL_N1}{Measurement of the right eye, Macula, Inner Nuclear Layer, Nasal Inner Parafovea}
#'     \item{R_MAC_INL_N2}{Measurement of the right eye, Macula, Inner Nuclear Layer, Nasal Outer Parafovea}
#'     \item{R_MAC_INL_I1}{Measurement of the right eye, Macula, Inner Nuclear Layer, Inferior Inner Parafovea}
#'     \item{R_MAC_INL_I2}{Measurement of the right eye, Macula, Inner Nuclear Layer, Inferior Outer Parafovea}
#'     \item{R_MAC_INL_T1}{Measurement of the right eye, Macula, Inner Nuclear Layer, Temporal Inner Parafovea}
#'     \item{R_MAC_INL_T2}{Measurement of the right eye, Macula, Inner Nuclear Layer, Temporal Outer Parafovea}
#'     \item{R_MAC_INL_C}{Measurement of the right eye, Macula, Inner Nuclear Layer, Center Fovea}
#'     \item{R_MAC_OPL_S1}{Measurement of the right eye, Macula, Outer Plexiform Layer, Superior Inner Parafovea}
#'     \item{R_MAC_OPL_S2}{Measurement of the right eye, Macula, Outer Plexiform Layer, Superior Outer Parafovea}
#'     \item{R_MAC_OPL_N1}{Measurement of the right eye, Macula, Outer Plexiform Layer, Nasal Inner Parafovea}
#'     \item{R_MAC_OPL_N2}{Measurement of the right eye, Macula, Outer Plexiform Layer, Nasal Outer Parafovea}
#'     \item{R_MAC_OPL_I1}{Measurement of the right eye, Macula, Outer Plexiform Layer, Inferior Inner Parafovea}
#'     \item{R_MAC_OPL_I2}{Measurement of the right eye, Macula, Outer Plexiform Layer, Inferior Outer Parafovea}
#'     \item{R_MAC_OPL_T1}{Measurement of the right eye, Macula, Outer Plexiform Layer, Temporal Inner Parafovea}
#'     \item{R_MAC_OPL_T2}{Measurement of the right eye, Macula, Outer Plexiform Layer, Temporal Outer Parafovea}
#'     \item{R_MAC_OPL_C}{Measurement of the right eye, Macula, Outer Plexiform Layer, Center Fovea}
#'     \item{R_MAC_ONL_S1}{Measurement of the right eye, Macula, Outer Nuclear Layer, Superior Inner Parafovea}
#'     \item{R_MAC_ONL_S2}{Measurement of the right eye, Macula, Outer Nuclear Layer, Superior Outer Parafovea}
#'     \item{R_MAC_ONL_N1}{Measurement of the right eye, Macula, Outer Nuclear Layer, Nasal Inner Parafovea}
#'     \item{R_MAC_ONL_N2}{Measurement of the right eye, Macula, Outer Nuclear Layer, Nasal Outer Parafovea}
#'     \item{R_MAC_ONL_I1}{Measurement of the right eye, Macula, Outer Nuclear Layer, Inferior Inner Parafovea}
#'     \item{R_MAC_ONL_I2}{Measurement of the right eye, Macula, Outer Nuclear Layer, Inferior Outer Parafovea}
#'     \item{R_MAC_ONL_T1}{Measurement of the right eye, Macula, Outer Nuclear Layer, Temporal Inner Parafovea}
#'     \item{R_MAC_ONL_T2}{Measurement of the right eye, Macula, Outer Nuclear Layer, Temporal Outer Parafovea}
#'     \item{R_MAC_ONL_C}{Measurement of the right eye, Macula, Outer Nuclear Layer, Center Fovea}
#'     \item{R_MAC_RPE_S1}{Measurement of the right eye, Macula, Retinal Pigment Epithelium, Superior Inner Parafovea}
#'     \item{R_MAC_RPE_S2}{Measurement of the right eye, Macula, Retinal Pigment Epithelium, Superior Outer Parafovea}
#'     \item{R_MAC_RPE_N1}{Measurement of the right eye, Macula, Retinal Pigment Epithelium, Nasal Inner Parafovea}
#'     \item{R_MAC_RPE_N2}{Measurement of the right eye, Macula, Retinal Pigment Epithelium, Nasal Outer Parafovea}
#'     \item{R_MAC_RPE_I1}{Measurement of the right eye, Macula, Retinal Pigment Epithelium, Inferior Inner Parafovea}
#'     \item{R_MAC_RPE_I2}{Measurement of the right eye, Macula, Retinal Pigment Epithelium, Inferior Outer Parafovea}
#'     \item{R_MAC_RPE_T1}{Measurement of the right eye, Macula, Retinal Pigment Epithelium, Temporal Inner Parafovea}
#'     \item{R_MAC_RPE_T2}{Measurement of the right eye, Macula, Retinal Pigment Epithelium, Temporal Outer Parafovea}
#'     \item{R_MAC_RPE_C}{Measurement of the right eye, Macula, Retinal Pigment Epithelium, Center Fovea}
#'     \item{R_MAC_PHOTO_S1}{Measurement of the right eye, Macula, Unknown Layer, Superior Inner Parafovea}
#'     \item{R_MAC_PHOTO_S2}{Measurement of the right eye, Macula, Unknown Layer, Superior Outer Parafovea}
#'     \item{R_MAC_PHOTO_N1}{Measurement of the right eye, Macula, Unknown Layer, Nasal Inner Parafovea}
#'     \item{R_MAC_PHOTO_N2}{Measurement of the right eye, Macula, Unknown Layer, Nasal Outer Parafovea}
#'     \item{R_MAC_PHOTO_I1}{Measurement of the right eye, Macula, Unknown Layer, Inferior Inner Parafovea}
#'     \item{R_MAC_PHOTO_I2}{Measurement of the right eye, Macula, Unknown Layer, Inferior Outer Parafovea}
#'     \item{R_MAC_PHOTO_T1}{Measurement of the right eye, Macula, Unknown Layer, Temporal Inner Parafovea}
#'     \item{R_MAC_PHOTO_T2}{Measurement of the right eye, Macula, Unknown Layer, Temporal Outer Parafovea}
#'     \item{R_MAC_PHOTO_C}{Measurement of the right eye, Macula, Unknown Layer, Center Fovea}
#'     \item{L_MAC_FULL_S1}{Measurement of the left eye, Macula, Complete Retinal Thickness, Superior Inner Parafovea}
#'     \item{L_MAC_FULL_S2}{Measurement of the left eye, Macula, Complete Retinal Thickness, Superior Outer Parafovea}
#'     \item{L_MAC_FULL_N1}{Measurement of the left eye, Macula, Complete Retinal Thickness, Nasal Inner Parafovea}
#'     \item{L_MAC_FULL_N2}{Measurement of the left eye, Macula, Complete Retinal Thickness, Nasal Outer Parafovea}
#'     \item{L_MAC_FULL_I1}{Measurement of the left eye, Macula, Complete Retinal Thickness, Inferior Inner Parafovea}
#'     \item{L_MAC_FULL_I2}{Measurement of the left eye, Macula, Complete Retinal Thickness, Inferior Outer Parafovea}
#'     \item{L_MAC_FULL_T1}{Measurement of the left eye, Macula, Complete Retinal Thickness, Temporal Inner Parafovea}
#'     \item{L_MAC_FULL_T2}{Measurement of the left eye, Macula, Complete Retinal Thickness, Temporal Outer Parafovea}
#'     \item{L_MAC_FULL_C}{Measurement of the left eye, Macula, Complete Retinal Thickness, Center Fovea}
#'     \item{L_MAC_RNFL_S1}{Measurement of the left eye, Macula, Retinal Nerve Fiber Layer, Superior Inner Parafovea}
#'     \item{L_MAC_RNFL_S2}{Measurement of the left eye, Macula, Retinal Nerve Fiber Layer, Superior Outer Parafovea}
#'     \item{L_MAC_RNFL_N1}{Measurement of the left eye, Macula, Retinal Nerve Fiber Layer, Nasal Inner Parafovea}
#'     \item{L_MAC_RNFL_N2}{Measurement of the left eye, Macula, Retinal Nerve Fiber Layer, Nasal Outer Parafovea}
#'     \item{L_MAC_RNFL_I1}{Measurement of the left eye, Macula, Retinal Nerve Fiber Layer, Inferior Inner Parafovea}
#'     \item{L_MAC_RNFL_I2}{Measurement of the left eye, Macula, Retinal Nerve Fiber Layer, Inferior Outer Parafovea}
#'     \item{L_MAC_RNFL_T1}{Measurement of the left eye, Macula, Retinal Nerve Fiber Layer, Temporal Inner Parafovea}
#'     \item{L_MAC_RNFL_T2}{Measurement of the left eye, Macula, Retinal Nerve Fiber Layer, Temporal Outer Parafovea}
#'     \item{L_MAC_RNFL_C}{Measurement of the left eye, Macula, Retinal Nerve Fiber Layer, Center Fovea}
#'     \item{L_MAC_GCL_S1}{Measurement of the left eye, Macula, Ganglion Cell Layer, Superior Inner Parafovea}
#'     \item{L_MAC_GCL_S2}{Measurement of the left eye, Macula, Ganglion Cell Layer, Superior Outer Parafovea}
#'     \item{L_MAC_GCL_N1}{Measurement of the left eye, Macula, Ganglion Cell Layer, Nasal Inner Parafovea}
#'     \item{L_MAC_GCL_N2}{Measurement of the left eye, Macula, Ganglion Cell Layer, Nasal Outer Parafovea}
#'     \item{L_MAC_GCL_I1}{Measurement of the left eye, Macula, Ganglion Cell Layer, Inferior Inner Parafovea}
#'     \item{L_MAC_GCL_I2}{Measurement of the left eye, Macula, Ganglion Cell Layer, Inferior Outer Parafovea}
#'     \item{L_MAC_GCL_T1}{Measurement of the left eye, Macula, Ganglion Cell Layer, Temporal Inner Parafovea}
#'     \item{L_MAC_GCL_T2}{Measurement of the left eye, Macula, Ganglion Cell Layer, Temporal Outer Parafovea}
#'     \item{L_MAC_GCL_C}{Measurement of the left eye, Macula, Ganglion Cell Layer, Center Fovea}
#'     \item{L_MAC_IPL_S1}{Measurement of the left eye, Macula, Inner Plexiform Layer, Superior Inner Parafovea}
#'     \item{L_MAC_IPL_S2}{Measurement of the left eye, Macula, Inner Plexiform Layer, Superior Outer Parafovea}
#'     \item{L_MAC_IPL_N1}{Measurement of the left eye, Macula, Inner Plexiform Layer, Nasal Inner Parafovea}
#'     \item{L_MAC_IPL_N2}{Measurement of the left eye, Macula, Inner Plexiform Layer, Nasal Outer Parafovea}
#'     \item{L_MAC_IPL_I1}{Measurement of the left eye, Macula, Inner Plexiform Layer, Inferior Inner Parafovea}
#'     \item{L_MAC_IPL_I2}{Measurement of the left eye, Macula, Inner Plexiform Layer, Inferior Outer Parafovea}
#'     \item{L_MAC_IPL_T1}{Measurement of the left eye, Macula, Inner Plexiform Layer, Temporal Inner Parafovea}
#'     \item{L_MAC_IPL_T2}{Measurement of the left eye, Macula, Inner Plexiform Layer, Temporal Outer Parafovea}
#'     \item{L_MAC_IPL_C}{Measurement of the left eye, Macula, Inner Plexiform Layer, Center Fovea}
#'     \item{L_MAC_INL_S1}{Measurement of the left eye, Macula, Inner Nuclear Layer, Superior Inner Parafovea}
#'     \item{L_MAC_INL_S2}{Measurement of the left eye, Macula, Inner Nuclear Layer, Superior Outer Parafovea}
#'     \item{L_MAC_INL_N1}{Measurement of the left eye, Macula, Inner Nuclear Layer, Nasal Inner Parafovea}
#'     \item{L_MAC_INL_N2}{Measurement of the left eye, Macula, Inner Nuclear Layer, Nasal Outer Parafovea}
#'     \item{L_MAC_INL_I1}{Measurement of the left eye, Macula, Inner Nuclear Layer, Inferior Inner Parafovea}
#'     \item{L_MAC_INL_I2}{Measurement of the left eye, Macula, Inner Nuclear Layer, Inferior Outer Parafovea}
#'     \item{L_MAC_INL_T1}{Measurement of the left eye, Macula, Inner Nuclear Layer, Temporal Inner Parafovea}
#'     \item{L_MAC_INL_T2}{Measurement of the left eye, Macula, Inner Nuclear Layer, Temporal Outer Parafovea}
#'     \item{L_MAC_INL_C}{Measurement of the left eye, Macula, Inner Nuclear Layer, Center Fovea}
#'     \item{L_MAC_OPL_S1}{Measurement of the left eye, Macula, Outer Plexiform Layer, Superior Inner Parafovea}
#'     \item{L_MAC_OPL_S2}{Measurement of the left eye, Macula, Outer Plexiform Layer, Superior Outer Parafovea}
#'     \item{L_MAC_OPL_N1}{Measurement of the left eye, Macula, Outer Plexiform Layer, Nasal Inner Parafovea}
#'     \item{L_MAC_OPL_N2}{Measurement of the left eye, Macula, Outer Plexiform Layer, Nasal Outer Parafovea}
#'     \item{L_MAC_OPL_I1}{Measurement of the left eye, Macula, Outer Plexiform Layer, Inferior Inner Parafovea}
#'     \item{L_MAC_OPL_I2}{Measurement of the left eye, Macula, Outer Plexiform Layer, Inferior Outer Parafovea}
#'     \item{L_MAC_OPL_T1}{Measurement of the left eye, Macula, Outer Plexiform Layer, Temporal Inner Parafovea}
#'     \item{L_MAC_OPL_T2}{Measurement of the left eye, Macula, Outer Plexiform Layer, Temporal Outer Parafovea}
#'     \item{L_MAC_OPL_C}{Measurement of the left eye, Macula, Outer Plexiform Layer, Center Fovea}
#'     \item{L_MAC_ONL_S1}{Measurement of the left eye, Macula, Outer Nuclear Layer, Superior Inner Parafovea}
#'     \item{L_MAC_ONL_S2}{Measurement of the left eye, Macula, Outer Nuclear Layer, Superior Outer Parafovea}
#'     \item{L_MAC_ONL_N1}{Measurement of the left eye, Macula, Outer Nuclear Layer, Nasal Inner Parafovea}
#'     \item{L_MAC_ONL_N2}{Measurement of the left eye, Macula, Outer Nuclear Layer, Nasal Outer Parafovea}
#'     \item{L_MAC_ONL_I1}{Measurement of the left eye, Macula, Outer Nuclear Layer, Inferior Inner Parafovea}
#'     \item{L_MAC_ONL_I2}{Measurement of the left eye, Macula, Outer Nuclear Layer, Inferior Outer Parafovea}
#'     \item{L_MAC_ONL_T1}{Measurement of the left eye, Macula, Outer Nuclear Layer, Temporal Inner Parafovea}
#'     \item{L_MAC_ONL_T2}{Measurement of the left eye, Macula, Outer Nuclear Layer, Temporal Outer Parafovea}
#'     \item{L_MAC_ONL_C}{Measurement of the left eye, Macula, Outer Nuclear Layer, Center Fovea}
#'     \item{L_MAC_RPE_S1}{Measurement of the left eye, Macula, Retinal Pigment Epithelium, Superior Inner Parafovea}
#'     \item{L_MAC_RPE_S2}{Measurement of the left eye, Macula, Retinal Pigment Epithelium, Superior Outer Parafovea}
#'     \item{L_MAC_RPE_N1}{Measurement of the left eye, Macula, Retinal Pigment Epithelium, Nasal Inner Parafovea}
#'     \item{L_MAC_RPE_N2}{Measurement of the left eye, Macula, Retinal Pigment Epithelium, Nasal Outer Parafovea}
#'     \item{L_MAC_RPE_I1}{Measurement of the left eye, Macula, Retinal Pigment Epithelium, Inferior Inner Parafovea}
#'     \item{L_MAC_RPE_I2}{Measurement of the left eye, Macula, Retinal Pigment Epithelium, Inferior Outer Parafovea}
#'     \item{L_MAC_RPE_T1}{Measurement of the left eye, Macula, Retinal Pigment Epithelium, Temporal Inner Parafovea}
#'     \item{L_MAC_RPE_T2}{Measurement of the left eye, Macula, Retinal Pigment Epithelium, Temporal Outer Parafovea}
#'     \item{L_MAC_RPE_C}{Measurement of the left eye, Macula, Retinal Pigment Epithelium, Center Fovea}
#'     \item{L_MAC_PHOTO_S1}{Measurement of the left eye, Macula, Unknown Layer, Superior Inner Parafovea}
#'     \item{L_MAC_PHOTO_S2}{Measurement of the left eye, Macula, Unknown Layer, Superior Outer Parafovea}
#'     \item{L_MAC_PHOTO_N1}{Measurement of the left eye, Macula, Unknown Layer, Nasal Inner Parafovea}
#'     \item{L_MAC_PHOTO_N2}{Measurement of the left eye, Macula, Unknown Layer, Nasal Outer Parafovea}
#'     \item{L_MAC_PHOTO_I1}{Measurement of the left eye, Macula, Unknown Layer, Inferior Inner Parafovea}
#'     \item{L_MAC_PHOTO_I2}{Measurement of the left eye, Macula, Unknown Layer, Inferior Outer Parafovea}
#'     \item{L_MAC_PHOTO_T1}{Measurement of the left eye, Macula, Unknown Layer, Temporal Inner Parafovea}
#'     \item{L_MAC_PHOTO_T2}{Measurement of the left eye, Macula, Unknown Layer, Temporal Outer Parafovea}
#'     \item{L_MAC_PHOTO_C}{Measurement of the left eye, Macula, Unknown Layer, Center Fovea}
#'     }
#'
#' Layer abbreviations include RNFL (Retinal Nerve Fiber Layer), GCL (Ganglion Cell Layer), IPL (Inner Plexiform Layer), INL (Inner Nuclear Layer), OPL (Outer Plexiform Layer), ONL (Outer Nuclear Layer), RPE (Retinal Pigment Epithelium), and IRL (Inner Retinal Layer).
#'
#' @source  Steingrube, N. (2023). Analysis of early changes in the retina and their association with diabetic alterations of the corneal nerve fiber plexus in type 2 diabetes mellitus. Unpublished doctoral dissertation. Faculty of Medicine, Heinrich-Heine University Dusseldorf.
#'
#' Department of Ophthalmology, University Hospital Dusseldorf, Heinrich Heine University, Germany
#'
#' @keywords datasets
#'
#'

NULL


#' Meta-analysis: Real World Evidence in metastatic colorectal cancer, comparing antiangiogenic treatments with chemotherapy
#'
#'
#' @name colon_cancer
#'
#' @docType data
#'
#' @description
#'
#' Meta-analysis of 7 RCTs, 4 cRWE studies, and 2 matched sRWE studies evaluating
#' progression-free survival (PFS) as a surrogate endpoint to overall survival (OS)
#' in metastatic colorectal cancer (mCRC), comparing antiangiogenic treatments with chemotherapy.
#'
#' @format  A dataframe with 13 rows and 6 columns. Each row represents study results, the columns are:
#'    \describe{
#'     \item{study}{Author and year.}
#'     \item{study_type}{randomized clinical trial or comparative/single-arm real-world-evidence}
#'     \item{pfs}{logarithm of hazard ratios of progression-free survival}
#'     \item{se_pfs}{standard error of pfs}
#'     \item{os}{logarithm of hazard ratios of overall survival}
#'     \item{se_os}{standard error of os}
#'     }
#'
#' @source  Wheaton L, Papanikos A, Thomas A, Bujkiewicz S. Using Bayesian Evidence Synthesis Methods to Incorporate Real-World Evidence in Surrogate Endpoint Evaluation. Medical Decision Making. 2023;43(5):539-552. doi:10.1177/0272989X231162852
#'
#' @keywords datasets
#'
#'

NULL



#' Meta-Analysis: Variation in False-Negative Rate of Reverse Transcriptase Polymerase Chain Reaction–Based SARS-CoV-2 Tests by Time Since Exposure
#'
#' @name fnrpcr
#'
#' @docType data
#'
#' @description
#' A dataset summarizing the variation in false-negative rates of reverse transcriptase polymerase chain reaction (RT-PCR)–based SARS-CoV-2 tests
#' as a function of time since exposure.
#'
#' @format A data frame with 410 rows and 11 columns. Each row represents the results from a study. The columns include:
#' \describe{
#'   \item{study}{Name of the author conducting the study.}
#'   \item{test}{Type of testing performed.}
#'   \item{day}{Number of days since symptom onset.}
#'   \item{day_min}{Minimum number of days since symptom onset. Applicable for studies by Guo et al. and Kim et al.}
#'   \item{day_max}{Maximum number of days since symptom onset. Applicable for studies by Guo et al. and Kim et al.}
#'   \item{n}{Total number of tests conducted on a given day.}
#'   \item{test_pos}{Number of positive test results.}
#'   \item{inconclusive}{Number of inconclusive test results. Applicable for studies by Kujawski et al. and Danis et al.}
#'   \item{nqp}{Number of positive but non-quantifiable test results, where the viral load is below the quantification threshold of log10(1) copies/1000 cells.}
#'   \item{pct_pos}{Proportion of positive tests expressed as a percentage.}
#' }
#'
#' @source Kucirka LM, Lauer SA, Laeyendecker O, Boon D, Lessler J. Variation in False-Negative Rate of Reverse Transcriptase Polymerase Chain Reaction-Based SARS-CoV-2 Tests by Time Since Exposure. Ann Intern Med. 2020 Aug 18;173(4):262-267. doi: 10.7326/M20-1495. Epub 2020 May 13. PMID: 32422057; PMCID: PMC7240870.
#'
#' @keywords datasets
#'
NULL



#' Meta-analysis: 83 observational studies assessing the effectiveness of intravitreal therapy for diabetic maculaedema
#'
#'
#' @name macula_rwe
#'
#' @docType data
#'
#' @description
#'
#' Meta-analysis of 83 studies comparing 12-month visual acuity change results in routine clinical practices
#' of intravitreal therapy for diabetic maculaedema (DME) to the change in RCTs by pooling data published in the
#' last decade on treated eyes of treatment effect from OS.
#'
#' @format  A dataframe with 83 rows and 13 columns. Each row represents study results, the columns are:
#'    \describe{
#'     \item{therapy}{Used Medication}
#'     \item{author_year}{Author and year.}
#'     \item{eyes}{Number of tested eyes.}
#'     \item{TE}{Mean Change in Visual Acuity after 12-Months. Where Visual Acuity is measured in VAR (Visual Acuity Rating) score}
#'     \item{seTE}{Standard Error of the Treatment Effect}
#'     \item{lower_95pct_ci}{Lower 95prc CI for TE}
#'     \item{upper_95pct_ci}{Upper 95prc CI for TE}
#'     \item{number_of_patients_at_baseline}{Number of Patients in Study at Baseline}
#'     \item{mean_age}{The mean age of patients per study}
#'     \item{baseline_va}{Mean Visual Acuity at Baseline. Where Visual Acuity is measured in VAR (Visual Acuity Rating) score}
#'     \item{12_month_va}{Mean Visual Acuity after 12 months. Where Visual Acuity is measured in VAR (Visual Acuity Rating) score}
#'     \item{baseline_cst}{Mean Central Subfield Thickness at Baseline}
#'     \item{12_month_cst}{Mean Central Subfield Thickness after 12 months}
#'     }
#'
#' @source Mehta H, Nguyen V, Barthelmes D, Pershing S, Chi GC, Dopart P, Gillies MC. Outcomes of Over 40,000 Eyes Treated for Diabetic Macula Edema in Routine Clinical Practice: A Systematic Review and Meta-analysis. Adv Ther. 2022 Dec;39(12):5376-5390. doi: 10.1007/s12325-022-02326-8. Epub 2022 Oct 15. PMID: 36241963; PMCID: PMC9618488.
#'
#'
#'
#' @keywords datasets
#'
#'

NULL

#' Meta-analysis: 29 randomized controlled studies (RCT) assessing the efficacy of acupuncture
#' treatments as complementary treatment in depression patients
#'
#'
#' @name acupuncture
#'
#' @docType data
#'
#' @description
#'
#' Meta-analysis of 29 studies on the effect of different methods of acupuncture
#' Therapy for depression compared to usual care control groups by pooling data from RCTs.
#'
#' @format  A dataframe with 29 rows and 11 columns. Each row represents study results, the columns are:
#'    \describe{
#'     \item{author_year}{Author and year.}
#'     \item{hedges_g}{changes in severity between intervention and control groups calculated using Hedges´g statistic}
#'     \item{std_err}{Standard Error of Hedges´g}
#'     \item{intervention}{treatment administered}
#'     \item{comparison}{control group treatment}
#'     \item{country}{origin country of the study}
#'     \item{sample_size}{total amount of patients per study}
#'     \item{number_treatments}{number of treatments received per study}
#'     \item{variation_acupuncture_points}{fixed: same acupuncture points used at each session; semi-fixed: some points pre-defined, some selected on the basis of the diagnosis/symptoms (location and amount); individualised: location and amount of points selected on basis of the diagnosis/symptoms}
#'     \item{number_acupuncture_points}{amount of acupuncture points for fixed-points-studies}
#'     \item{NICMAN}{NICMAN scale Points to evaluate the Quality of the administered acupuncture}
#'     \item{random_sequence_generation}{Risk of selection bias (Random sequence generation) low risk of bias: high, high risk: low, unclear: unclear}
#'     \item{allocation_concealment}{Risk of selection bias (allocation concealment) low risk of bias: high, high risk: low, unclear: unclear}
#'     \item{blinding_participants_personnel}{Risk of performance bias (blinding of participants and personnel) low risk of bias: high, high risk: low, unclear: unclear}
#'     \item{blinding_outcome_assessment}{Risk of detection bias (blinding oft outcome assessment) low risk of bias: high, high risk: low, unclear: unclear}
#'     \item{incomplete_outcome_data}{Risk of attrition bias (incomplete outcome data) low risk of bias: high, high risk: low, unclear: unclear}
#'     \item{selective_reporting}{Risk of reporting bias (selective reporting) low risk of bias: high, high risk: low, unclear: unclear}
#'     \item{other_bias}{Risk of other biases; low risk of bias: high, high risk: low, unclear: unclear}
#'     }
#'
#'
#' @source  Armour M, Smith CA, Wang LQ, Naidoo D, Yang GY, MacPherson H, Lee MS, Hay P. Acupuncture for Depression: A Systematic Review and Meta-Analysis. J Clin Med. 2019 Jul 31;8(8):1140. doi: 10.3390/jcm8081140. PMID: 31370200; PMCID: PMC6722678.
#'
#'
#'
#' @keywords datasets
#'
#'

NULL


#' Meta-analysis: generalized evidence synthesis of total hip replacement
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
#' Meta-analysis of 35 Observational Studies from PubMed, Cocharane Library and SciELO databases
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

#' Meta-analysis: 31 randomized controled trials (RCTs) with reported discrepancies
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
#' @format  A matrix with 31 rows and 11 columns. Each row represents study results, the columns are:
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
#'     \item{Sequence}{Bias arising from the randomization process.}
#'     \item{Allocation}{Bias due to deviations from intended interventions.}
#'     \item{Blinding}{Bias introduced by lack of blinding.}
#'     \item{Outcome}{Bias in measurement of the outcome.}
#'     \item{Reporting}{Bias in selection of the reported result.}
#'     \item{Other}{Selection bias, performance bias, detection bias, attrition bias, etc.}
#'
#'     }
#'
#' @references Verde, P. E. (2017) Two Examples of Bayesian Evidence Synthesis with the Hierarchical Meta-Regression Approach. Chap.9, pag 189-206. Bayesian Inference, ed. Tejedor, Javier Prieto. InTech.
#'
#'
#' @source  Nowbar, A N, et al. (2014) Discrepancies in autologous bone marrow stem cell trials and enhancement of ejection fraction (DAMASCENE): weighted regression and meta-analysis. BMJ, 348,1-9.
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

