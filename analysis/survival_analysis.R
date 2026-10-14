SHAP_ONLY <- FALSE

required_packages <- c(
  "broom",
  "caret",
  "dplyr",
  "forcats",
  "future",
  "ggDCA",
  "ggExtra",
  "ggplot2",
  "ggprism",
  "ggpubr",
  "ggsignif",
  "gridExtra",
  "haven",
  "Hmisc",
  "lubridate",
  "MASS",
  "mlr3extralearners",
  "mlr3proba",
  "mlr3verse",
  "nomogramFormula",
  "RColorBrewer",
  "readxl",
  "recipes",
  "reshape",
  "reshape2",
  "riskRegression",
  "rms",
  "rsample",
  "stringr",
  "survIDINRI",
  "survival",
  "survivalsvm",
  "survminer",
  "tableone",
  "tibble",
  "tidyr",
  "timeROC",
  "tidyverse"
)
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop("Install required packages: ", paste(missing_packages, collapse = ", "))
}
invisible(lapply(required_packages, library, character.only = TRUE))

project_dir <- normalizePath(
  Sys.getenv("SAGITTAL_FD_PROJECT_DIR", unset = "."),
  mustWork = FALSE
)
data_dir <- Sys.getenv(
  "SAGITTAL_FD_DATA_DIR",
  unset = file.path(project_dir, "data")
)
output_dir <- Sys.getenv(
  "SAGITTAL_FD_OUTPUT_DIR",
  unset = file.path(project_dir, "results")
)
development_clinical_file <- Sys.getenv(
  "SAGITTAL_FD_DEVELOPMENT_CLINICAL",
  unset = file.path(data_dir, "development_clinical.sav")
)
development_fd_file <- Sys.getenv(
  "SAGITTAL_FD_DEVELOPMENT_FEATURES",
  unset = file.path(data_dir, "development_fd_features.csv")
)
external_clinical_file <- Sys.getenv(
  "SAGITTAL_FD_EXTERNAL_CLINICAL",
  unset = file.path(data_dir, "external_clinical.xlsx")
)
external_fd_file <- Sys.getenv(
  "SAGITTAL_FD_EXTERNAL_FEATURES",
  unset = file.path(data_dir, "external_fd_features.csv")
)

input_files <- c(
  development_clinical_file,
  development_fd_file,
  external_clinical_file,
  external_fd_file
)
missing_inputs <- input_files[!file.exists(input_files)]
if (length(missing_inputs)) {
  stop(
    "Missing input files. Configure the SAGITTAL_FD_* environment variables: ",
    paste(missing_inputs, collapse = ", ")
  )
}

model_dir <- file.path(output_dir, "models")
figure_dir <- file.path(output_dir, "figures")
explain_dir <- file.path(output_dir, "explainability")
dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(explain_dir, recursive = TRUE, showWarnings = FALSE)

extract_survival_probabilities <- function(prediction, timepoints, n_observations) {
  # Keep the full distribution: subsetting may drop statistical decorators.
  distributions <- prediction$distr
  survival_function <- distributions$survival
  if (!is.function(survival_function)) {
    cdf_function <- distributions$cdf
    if (!is.function(cdf_function)) {
      stop("Prediction distribution has neither survival() nor cdf(): ",
           paste(class(distributions), collapse = "/"))
    }
    survival_function <- function(time) 1 - cdf_function(time)
  }
  if (!is.numeric(timepoints) || !length(timepoints) ||
      any(!is.finite(timepoints)) || any(timepoints < 0)) {
    stop("timepoints must contain finite, non-negative numeric times.")
  }
  # One time per call makes orientation unambiguous, even for square matrices.
  probabilities <- vapply(timepoints, function(time) {
    values <- as.numeric(survival_function(time))
    if (length(values) != n_observations) {
      stop("Expected ", n_observations, " patient probabilities; received ",
           length(values), " at time ", time, ".")
    }
    if (any(!is.finite(values)) || any(values < 0 | values > 1)) {
      stop("Invalid survival probabilities at time ", time, ".")
    }
    values
  }, numeric(n_observations))
  probabilities <- matrix(probabilities, nrow = n_observations,
                          ncol = length(timepoints))
  colnames(probabilities) <- as.character(timepoints)
  as.data.frame(probabilities, check.names = FALSE)
}

if (!SHAP_ONLY) {
source(file.path(project_dir, "R", "evaluation_helpers.R"))
itps <- c(6 * c(1, 2,3, 4,5,6,7,8,9,10,11,12))
Clinical_Data <- read_sav(development_clinical_file)
FD_Data <- read.csv(development_fd_file)
FD_Data$PatientsNum <- as.numeric(str_sub(FD_Data$patient_id,start = 0, end = 7))
merged_data <- merge(Clinical_Data, FD_Data[,-which(colnames(FD_Data) == "patient_id")], by = "PatientsNum")
merged_data$BM <- ifelse(merged_data$BM == 1, 1, 0)
merged_data$status_os <- ifelse(merged_data$status_os == 1, 1, 0)
merged_data[,c(2,4,12,15,39,41:42,45)] <- lapply(merged_data[,c(2,4,12,15,39,41:42,45)],function(x){
  factor(x,levels = c(data.frame(Hmisc::valueLabel(x))[,1]),labels = rownames(data.frame(Hmisc::valueLabel(x))))
})
non_pcidata <- merged_data
non_pcidata$ECOG <- ifelse(non_pcidata$ECOG>1,">1","0-1")
non_pcidata$PCI <- ifelse(non_pcidata$PCI==0,"No", "Yes")
non_pcidata$RT <- ifelse(non_pcidata$TimeToRT==0,"No", "Yes")
non_pcidata$Chemotherapy <- ifelse(non_pcidata$Chemotherapy>=4,"Yes", "No")

non_pcidata$PatientsNum_chr <- as.character(non_pcidata$PatientsNum)
non_pcidata$year <- as.integer(substr(non_pcidata$PatientsNum_chr, start = 1, stop = 2))+2000
biomarkers <- c(
  "LDH",  "NSE", "Cyfra211", "ProGRP"
)
MoLan_Data1 <- non_pcidata[c(1,2:4,8,15,40,50,39,67:82,64,63,85)]
MoLan_Data1 <- MoLan_Data1[!MoLan_Data1$TreatmentMode %in% c("CT_IM" ,"CT_IM_RT"),]
MoLan_Data1$TreatmentMode <- ifelse(MoLan_Data1$TreatmentMode %in% c("CRT","C_CRT_C"),"CRT",
                                    ifelse(MoLan_Data1$TreatmentMode=="C","C","CTRT"))
MoLan_Data1$group <- "shanxi"


FD_Data_tianjin <- read.csv(external_fd_file)
FD_Data_tianjin$PatientsNum <- as.numeric(str_extract(FD_Data_tianjin$patient_id, "^\\d+"))
clinical_data_tianjin <- read_xlsx(external_clinical_file)
clinical_data_tianjin$PatientsNum <- clinical_data_tianjin$序号
clinical_data_tianjin$ECOG <- ifelse(clinical_data_tianjin$`KPS≥80`<=1,"0-1",">1")
clinical_data_tianjin$Chemotherapy <- ifelse(clinical_data_tianjin$化疗总周期数<4,"No","Yes")
clinical_data_tianjin$RT <- "Yes"
clinical_data_tianjin$PCI <- ifelse(clinical_data_tianjin$PCI==1,"Yes","No")
clinical_data_tianjin <- merge.data.frame(clinical_data_tianjin,FD_Data_tianjin,by = "PatientsNum")

clinical_data_tianjin$status_os <- clinical_data_tianjin$状态
clinical_data_tianjin <- clinical_data_tianjin[!(clinical_data_tianjin$OS<6&clinical_data_tianjin$status_os==0),]
clinical_data_tianjin$Age <- clinical_data_tianjin$年龄
clinical_data_tianjin$Gender <- ifelse(clinical_data_tianjin$性别==1,"Male","Female")
clinical_data_tianjin$Clinical_Stage <- ifelse(clinical_data_tianjin$`临床分期（1:I期;2:II期;3:IIIA期;4:IIIB期;5:IIIC期）`==3,"IIIa",
                                               ifelse(clinical_data_tianjin$`临床分期（1:I期;2:II期;3:IIIA期;4:IIIB期;5:IIIC期）`==4,"IIIb",
                                                      ifelse(clinical_data_tianjin$`临床分期（1:I期;2:II期;3:IIIA期;4:IIIB期;5:IIIC期）`==5,"III‎c","I-IIb")))
clinical_data_tianjin$Smoke <- ifelse(clinical_data_tianjin$吸烟==1,"Yes","No")
clinical_data_tianjin$TreatmentMode <- ifelse(clinical_data_tianjin$同步化疗有无==1, "CRT","CTRT")
clinical_data_tianjin$ITHScore <- clinical_data_tianjin$ith_score
clinical_data_tianjin$date_correct_base <- as.Date(clinical_data_tianjin$确诊时间, format = "%m/%d/%Y")
clinical_data_tianjin$year <- as.integer(format(clinical_data_tianjin$date_correct_base, "%Y"))
clinical_data_tianjin$group <- "tianjin"

clinical_data_tianjin <- clinical_data_tianjin %>%
  mutate(
    NEU      = 治疗前中性粒细胞,
    LYM      = 治疗前淋巴细胞,
    HGB      = 治疗前血红蛋白,
    Alb      = 治疗前ALB,
    LDH      = 治疗前LDH,
    CEA      = 治疗前CEA,
    CA199    = 治疗前CA199,
    NSE      = 治疗前NSE,
    Cyfra211 = 治疗前Cyfra21,
    ProGRP   = 治疗前ProGRP
  )

MoLan_Data2 <- clinical_data_tianjin[names(MoLan_Data1)]
MoLan_Data <- rbind.data.frame(MoLan_Data1,MoLan_Data2)
MoLan_Data$group <- ifelse(MoLan_Data$year>2018,">2018","2013-2018")
MoLan_Data <- na.omit(MoLan_Data)
seednum <- 377
set.seed(377)
datasplit <- rsample::initial_split(MoLan_Data, prop = 0.8, strata = "status_os")
temp_data <- rsample::training(datasplit)
validation_data2 <- rsample::testing(datasplit)
split_2 <- rsample::initial_split(temp_data, prop = 0.7, strata = "status_os")
train_data <- rsample::training(split_2)
validation_data1 <- rsample::testing(split_2)

train_data$group <- "train"
validation_data1$group <- "validation1"
validation_data2$group <- "validation2"
MoLan_Data2 <- bind_rows(train_data,validation_data1,validation_data2)
data_Unweighted <- CreateTableOne(vars = names(MoLan_Data2)[c(2:9,29)],addOverall = T,factorVars =  names(MoLan_Data2)[c(2,4:9,29)],strata = "group",
                                  argsApprox = F,
                                  data = MoLan_Data2)
tab3mat <- print(data_Unweighted,showAllLevels = T,smd = T)
mydata <- data.frame(tab3mat)
survfit(Surv(OS,status_os==0)~1,data = MoLan_Data2)
UnicoX <- function(x, data){
  FML <- as.formula(paste0("Surv(OS, status_os) ~ ", x))
  coxph_model <- survival::coxph(FML, data = data)
  diff_data <- summary(coxph_model)
  HR <- round(diff_data$coefficients[,2], 3)
  HR_L <- round(diff_data$conf.int[,3], 3)
  HR_U <- round(diff_data$conf.int[,4], 3)
  Unicox <- data.frame(
    "Characteristics" = rownames(diff_data$coefficients),
    "HR(95%CI)" = paste(HR, "(", HR_L, "-", HR_U, ")"),
    "P value" = round(diff_data$coefficients[,"Pr(>|z|)"], 3),
    "var_name" = x
  )
  return(Unicox)
}
uni_names <- colnames(train_data)[c(2:25)]
uni_univar <- lapply(uni_names, UnicoX, data = train_data)
uni_univar <- do.call(rbind, uni_univar)

sig_vars <- unique(uni_univar[uni_univar$P.value<0.05,"var_name"])

dd <- datadist(train_data)
options(datadist = "dd")
FML <- as.formula(paste0("Surv(OS, status_os) ~ ", paste(unique(sig_vars), collapse = " + ")))
MultiCox <- stats::step(coxph(FML,data = train_data),trace = 0)
diff_data <- summary(MultiCox)
MultiCoxname <- rownames(data.frame(diff_data$coefficients))
Pvalue=round(diff_data$coefficients[,5],3)
HR=round(diff_data$coefficients[,2],3)
HR_L=round(diff_data$conf.int[,3],3)
HR_H=round(diff_data$conf.int[,4],3)
HR_result=paste0(HR,"(",HR_L,"-",HR_H,")")
Mulcox <- data.frame("Characteristics"=MultiCoxname,
                     "Hazard Ratio (95%CI)"=HR_result,
                     "P value"=Pvalue)
Final_PFS_PSM <- merge.data.frame(uni_univar,Mulcox,by="Characteristics",all = T,sort = T)
Final_PFS_PSM
sig_vars <- unique(Final_PFS_PSM[!is.na(Final_PFS_PSM$P.value.y),"var_name"])


library(survival)
library(ggplot2)

# 提取已拟合 Cox 模型的系数和协方差矩阵
beta <- coef(MultiCox)
stopifnot(all(is.finite(beta)))

V <- vcov(MultiCox)
V <- V[names(beta), names(beta), drop = FALSE]
R <- cov2cor(V)

# 将各系数对应回原始变量
X <- model.matrix(MultiCox)
term_names <- attr(terms(MultiCox), "term.labels")
column_index <- match(names(beta), colnames(X))

stopifnot(!anyNA(column_index))
term_index <- attr(X, "assign")[column_index]
stopifnot(length(term_index) == length(beta))

# 正定性检查，并计算完整矩阵的对数行列式
logdet_R <- 2 * sum(log(diag(chol(R))))

vif_data <- data.frame(
  Variable = character(),
  Df = integer(),
  GVIF = numeric(),
  VIF = numeric()
)

# 按原始变量计算，分类变量的多个系数作为一组
for (k in sort(unique(term_index))) {

  selected <- which(term_index == k)
  remaining <- setdiff(seq_along(beta), selected)
  df <- length(selected)

  R_selected <- R[selected, selected, drop = FALSE]
  logdet_selected <- 2 * sum(log(diag(chol(R_selected))))

  logdet_remaining <- 0

  if (length(remaining) > 0) {
    R_remaining <- R[remaining, remaining, drop = FALSE]
    logdet_remaining <- 2 * sum(log(diag(chol(R_remaining))))
  }

  gvif <- exp(logdet_selected + logdet_remaining - logdet_R)

  vif_data <- rbind(
    vif_data,
    data.frame(
      Variable = term_names[k],
      Df = df,
      GVIF = gvif,
      VIF = gvif^(1 / df)
    )
  )
}

# 同时保留标准的调整后 GVIF，方便报告
vif_data$Adjusted_GVIF <- vif_data$GVIF^(1 / (2 * vif_data$Df))

vif_data <- vif_data[order(vif_data$VIF, decreasing = TRUE), ]
rownames(vif_data) <- NULL

print(vif_data)

vif_data$Variable <- factor(
  vif_data$Variable,
  levels = vif_data$Variable
)

# 默认上限为10；若超过10，自动扩大，避免隐藏柱子
y_max <- max(10, ceiling(max(vif_data$VIF) * 1.15))

ggplot(vif_data, aes(x = Variable, y = VIF)) +
  geom_col(fill = "#00468B", width = 0.7) +
  geom_hline(
    yintercept = 5,
    color = "#E60049",
    linetype = "dashed",
    linewidth = 1
  ) +
  geom_text(
    aes(label = round(VIF, 2)),
    vjust = -0.3,
    size = 3.5,
    color = "black"
  ) +
  theme_bw() +
  theme(
    text = element_text(
      family = ifelse(
        Sys.info()["sysname"] == "Windows",
        "SimHei", "Arial"
      ),
      size = 10
    ),
    axis.title = element_text(face = "bold", size = 11),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y = element_text(size = 9),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    panel.border = element_rect(color = "black", linewidth = 0.8),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12)
  ) +
  labs(
    x = "Variables",
    y = expression(GVIF^{1/Df}),
    title = "GVIF Analysis for Multicollinearity"
  ) +
  coord_cartesian(ylim = c(0, y_max))


traindata <- train_data[,c(sig_vars,"OS","status_os")]
validation_data2_model_data <- validation_data2[,c(sig_vars,"OS","status_os")]
validation_data1 <- validation_data1[,c(sig_vars,"OS","status_os")]

cat("\n========== 开始训练 Cox 模型 ==========\n")

datarecipe_coxph <- recipe(OS + status_os ~ ., traindata) %>%
  prep()
datarecipe_coxph

traindata2_coxph <- bake(datarecipe_coxph, new_data = NULL) %>%
  dplyr::select(OS, status_os, everything())
validation_data2_coxph <- bake(datarecipe_coxph, new_data = validation_data2_model_data) %>%
  dplyr::select(OS, status_os, everything())

task_train_coxph <- as_task_surv(
  traindata2_coxph,
  time = "OS",
  event = "status_os",
  type = "right"
)
task_validation2_coxph <- as_task_surv(
  validation_data2_coxph,
  time = "OS",
  event = "status_os",
  type = "right"
)

learner_coxph <- lrn("surv.coxph")
learner_coxph

set.seed(seednum)
learner_coxph$train(task_train_coxph)
learner_coxph$model
summary(learner_coxph$model)

predtrain_coxph <- learner_coxph$predict(task_train_coxph)
predtrain_coxph$score(msrs(c("surv.cindex")))

predprobtrain_coxph <-
  predtrain_coxph$distr[1:nrow(traindata2_coxph)]$survival(itps) %>%
  t() %>%
  as.data.frame() %>%
  mutate(model = "coxph",
         dataset = "train",
         time = traindata2_coxph$OS,
         status_os = traindata2_coxph$status_os)

evaltrain_coxph <- eval4sa(
  predprob = predprobtrain_coxph,
  preddata = traindata2_coxph,
  etime = "OS",
  estatus = "status_os",
  model = "coxph",
  dataset = "train",
  timepoints = itps,
  plotcalimethod = "quantile",
  bw4nne = NULL,
  q4quantile = 1
)

predvalidation2_coxph <- learner_coxph$predict(task_validation2_coxph)
predprobvalidation2_coxph <-
  predvalidation2_coxph$distr[1:nrow(validation_data2_coxph)]$survival(itps) %>%
  t() %>%
  as.data.frame() %>%
  mutate(model = "coxph",
         dataset = "validation2",
         time = validation_data2_coxph$OS,
         status_os = validation_data2_coxph$status_os)

evalvalidation2_coxph <- eval4sa(
  predprob = predprobvalidation2_coxph,
  preddata = validation_data2_coxph,
  etime = "OS",
  estatus = "status_os",
  model = "coxph",
  dataset = "validation2",
  timepoints = itps,
  plotcalimethod = "quantile",
  bw4nne = NULL,
  q4quantile = 1
)
evaltrain_coxph$auc
evalvalidation2_coxph$auc
save(predtrain_coxph, predprobtrain_coxph, evaltrain_coxph,
     predvalidation2_coxph, predprobvalidation2_coxph, evalvalidation2_coxph,
     file = file.path(model_dir, "coxph.RData"))

cat("========== Cox 模型训练完成 ==========\n\n")

cat("\n========== 开始训练 Lasso 模型 ==========\n")

datarecipe_lasso <- recipe(OS + status_os ~ ., traindata) %>%
  step_dummy(all_nominal_predictors()) %>%
  prep()

traindata2_lasso <- bake(datarecipe_lasso, new_data = NULL) %>%
  dplyr::select(OS, status_os, everything())
validation_data2_lasso <- bake(datarecipe_lasso, new_data = validation_data2_model_data) %>%
  dplyr::select(OS, status_os, everything())

task_train_lasso <- as_task_surv(
  traindata2_lasso,
  time = "OS",
  event = "status_os",
  type = "right"
)
task_validation2_lasso <- as_task_surv(
  validation_data2_lasso,
  time = "OS",
  event = "status_os",
  type = "right"
)

learner_lasso_0 <- lrn(
  "surv.cv_glmnet",
  type.measure = "C",
  s = "lambda.1se"
)
learner_lasso_1 <- ppl(
  "distrcompositor",
  learner = learner_lasso_0,
  estimator = "kaplan",
  form = "ph"
)
learner_lasso <- as_learner(learner_lasso_1)

set.seed(seednum)
learner_lasso$train(task_train_lasso)
coef(learner_lasso$model$surv.cv_glmnet$model)

predtrain_lasso <- learner_lasso$predict(task_train_lasso)
predprobtrain_lasso <-
  predtrain_lasso$distr[1:nrow(traindata2_lasso)]$survival(itps) %>%
  t() %>%
  as.data.frame() %>%
  mutate(model = "lasso",
         dataset = "train",
         time = traindata2_lasso$OS,
         status = traindata2_lasso$status_os)

evaltrain_lasso <- eval4sa(
  predprob = predprobtrain_lasso,
  preddata = traindata2_lasso,
  etime = "OS",
  estatus = "status_os",
  model = "lasso",
  dataset = "train",
  timepoints = itps,
  plotcalimethod = "quantile",
  bw4nne = NULL,
  q4quantile = 1
)

predvalidation2_lasso <- learner_lasso$predict(task_validation2_lasso)
predprobvalidation2_lasso <-
  predvalidation2_lasso$distr[1:nrow(validation_data2_lasso)]$survival(itps) %>%
  t() %>%
  as.data.frame() %>%
  mutate(model = "lasso",
         dataset = "validation2",
         time = validation_data2_lasso$OS,
         status = validation_data2_lasso$status_os)

evalvalidation2_lasso <- eval4sa(
  predprob = predprobvalidation2_lasso,
  preddata = validation_data2_lasso,
  etime = "OS",
  estatus = "status_os",
  model = "lasso",
  dataset = "validation2",
  timepoints = itps,
  plotcalimethod = "quantile",
  bw4nne = NULL,
  q4quantile = 1
)
evaltrain_lasso$auc
evalvalidation2_lasso$auc
save(predtrain_lasso, predprobtrain_lasso, evaltrain_lasso,
     predvalidation2_lasso, predprobvalidation2_lasso, evalvalidation2_lasso,
     file = file.path(model_dir, "lasso.RData"))

cat("========== Lasso 模型训练完成 ==========\n\n")

cat("\n========== 开始训练 决策树 模型 ==========\n")

datarecipe_rpart <- recipe(OS + status_os ~ ., traindata) %>%
  prep()

traindata2_rpart <- bake(datarecipe_rpart, new_data = NULL) %>%
  dplyr::select(OS, status_os, everything())
validation_data2_rpart <- bake(datarecipe_rpart, new_data = validation_data2_model_data) %>%
  dplyr::select(OS, status_os, everything())

task_train_rpart <- as_task_surv(
  traindata2_rpart,
  time = "OS",
  event = "status_os",
  type = "right"
)
task_validation2_rpart <- as_task_surv(
  validation_data2_rpart,
  time = "OS",
  event = "status_os",
  type = "right"
)

learner_rpart_0 <- ppl(
  "distrcompositor",
  learner = lrn(
    "surv.rpart",
    cp = to_tune(0.001, 0.5),
    minbucket = to_tune(5, 9)
  ),
  estimator = "kaplan",
  form = "ph"
) %>%
  as_learner()

learner_rpart <- auto_tuner(
  tuner = tnr(
    "grid_search",
    param_resolutions = c(surv.rpart.cp = 5,
                          surv.rpart.minbucket = 3),
    batch_size = 4
  ),
  learner = learner_rpart_0,
  resampling = rsmp("cv", folds = 5),
  measure = msr("surv.cindex"),
  terminator = trm("none")
)

future::plan(future::sequential)

set.seed(seednum)
learner_rpart$train(task_train_rpart)
learner_rpart$tuning_result

learner_rpart$predict_type <- "crank"
baseline_rpart <- mlr3::lrn("surv.kaplan")
baseline_rpart$train(task_train_rpart)

compose_rpart_prediction <- function(task) {
  pred <- learner_rpart$predict(task)
  risk <- as.numeric(pred$crank)
  if (any(!is.finite(risk)) || any(risk <= 0)) {
    stop("rpart relative risks must be finite and positive for log transformation.")
  }
  pred_lp <- mlr3proba::PredictionSurv$new(
    row_ids = pred$row_ids,
    truth = pred$truth,
    crank = risk,
    lp = log(risk)
  )
  # Compose probabilities under PH using a training-only Kaplan baseline.
  compositor <- mlr3pipelines::po(
    "distrcompose",
    param_vals = list(form = "ph", overwrite = TRUE)
  )
  compositor$train(list(NULL, NULL))
  result <- compositor$predict(list(
    base = baseline_rpart$predict(task),
    pred = pred_lp
  ))[[1]]
  if (is.null(result$distr)) stop("rpart distribution composition failed.")
  result
}

predtrain_rpart <- compose_rpart_prediction(task_train_rpart)
predprobtrain_rpart <-
  extract_survival_probabilities(
    prediction = predtrain_rpart,
    timepoints = itps,
    n_observations = nrow(traindata2_rpart)
  ) %>%
  mutate(model = "rpart",
         dataset = "train",
         time = traindata2_rpart$OS,
         status = traindata2_rpart$status_os)

evaltrain_rpart <- eval4sa(
  predprob = predprobtrain_rpart,
  preddata = traindata2_rpart,
  etime = "OS",
  estatus = "status_os",
  model = "rpart",
  dataset = "train",
  timepoints = itps,
  plotcalimethod = "quantile",
  bw4nne = NULL,
  q4quantile = 1
)

predvalidation2_rpart <- compose_rpart_prediction(task_validation2_rpart)
predprobvalidation2_rpart <-
  extract_survival_probabilities(
    prediction = predvalidation2_rpart,
    timepoints = itps,
    n_observations = nrow(validation_data2_rpart)
  ) %>%
  mutate(model = "rpart",
         dataset = "validation2",
         time = validation_data2_rpart$OS,
         status = validation_data2_rpart$status_os)

evalvalidation2_rpart <- eval4sa(
  predprob = predprobvalidation2_rpart,
  preddata = validation_data2_rpart,
  etime = "OS",
  estatus = "status_os",
  model = "rpart",
  dataset = "validation2",
  timepoints = itps,
  plotcalimethod = "quantile",
  bw4nne = NULL,
  q4quantile = 1
)
evaltrain_rpart$auc
evalvalidation2_rpart$auc
save(predtrain_rpart, predprobtrain_rpart, evaltrain_rpart,
     predvalidation2_rpart, predprobvalidation2_rpart, evalvalidation2_rpart,
     file = file.path(model_dir, "ctree.RData"))

cat("========== 决策树 模型训练完成 ==========\n\n")

cat("\n========== 开始训练 随机森林 模型 ==========\n")

datarecipe_rsf <- recipe(OS + status_os ~ ., traindata) %>%
  prep()

traindata2_rsf <- bake(datarecipe_rsf, new_data = NULL) %>%
  dplyr::select(OS, status_os, everything())
validation_data2_rsf <- bake(datarecipe_rsf, new_data = validation_data2_model_data) %>%
  dplyr::select(OS, status_os, everything())

task_train_rsf <- as_task_surv(
  traindata2_rsf,
  time = "OS",
  event = "status_os",
  type = "right"
)
task_validation2_rsf <- as_task_surv(
  validation_data2_rsf,
  time = "OS",
  event = "status_os",
  type = "right"
)

learner_rsf_0 <- lrn(
  "surv.rfsrc",
  ntree = to_tune(200, 500),
  mtry = to_tune(3, 5),
  nodesize = to_tune(15, 21),
  predict_type = "distr"
)

learner_rsf <- auto_tuner(
  tuner = tnr("grid_search", resolution = 3, batch_size = 3),
  learner = learner_rsf_0,
  resampling = rsmp("cv", folds = 5),
  measure = msr("surv.cindex"),
  terminator = trm("none")
)

future::plan(future::sequential)

set.seed(seednum)
learner_rsf$train(task_train_rsf)
learner_rsf$tuning_result

predtrain_rsf <- learner_rsf$predict(task_train_rsf)
predprobtrain_rsf <-
  predtrain_rsf$distr[1:nrow(traindata2_rsf)]$survival(itps) %>%
  t() %>%
  as.data.frame() %>%
  mutate(model = "rsf",
         dataset = "train",
         time = traindata2_rsf$OS,
         status = traindata2_rsf$status_os)

evaltrain_rsf <- eval4sa(
  predprob = predprobtrain_rsf,
  preddata = traindata2_rsf,
  etime = "OS",
  estatus = "status_os",
  model = "rsf",
  dataset = "train",
  timepoints = itps,
  plotcalimethod = "quantile",
  bw4nne = NULL,
  q4quantile = 1
)

predvalidation2_rsf <- learner_rsf$predict(task_validation2_rsf)
predprobvalidation2_rsf <-
  predvalidation2_rsf$distr[1:nrow(validation_data2_rsf)]$survival(itps) %>%
  t() %>%
  as.data.frame() %>%
  mutate(model = "rsf",
         dataset = "validation2",
         time = validation_data2_rsf$OS,
         status = validation_data2_rsf$status_os)

evalvalidation2_rsf <- eval4sa(
  predprob = predprobvalidation2_rsf,
  preddata = validation_data2_rsf,
  etime = "OS",
  estatus = "status_os",
  model = "rsf",
  dataset = "validation2",
  timepoints = itps,
  plotcalimethod = "quantile",
  bw4nne = NULL,
  q4quantile = 1
)
save(predtrain_rsf, predprobtrain_rsf, evaltrain_rsf,
     predvalidation2_rsf, predprobvalidation2_rsf, evalvalidation2_rsf,
     file = file.path(model_dir, "rsf.RData"))

cat("========== 随机森林 模型训练完成 ==========\n\n")

cat("\n========== 开始训练 GBM 模型 ==========\n")

datarecipe_gbm <- recipe(OS + status_os ~ ., traindata) %>%
  step_dummy(all_nominal_predictors()) %>%
  prep()

traindata2_gbm <- bake(datarecipe_gbm, new_data = NULL) %>%
  dplyr::select(OS, status_os, everything())
validation_data2_gbm <- bake(datarecipe_gbm, new_data = validation_data2_model_data) %>%
  dplyr::select(OS, status_os, everything())

task_train_gbm <- as_task_surv(
  traindata2_gbm,
  time = "OS",
  event = "status_os",
  type = "right"
)
task_validation2_gbm <- as_task_surv(
  validation_data2_gbm,
  time = "OS",
  event = "status_os",
  type = "right"
)

learner_gbm_0 <- lrn(
  "surv.gbm",
  n.trees = to_tune(100L, 300L),
  interaction.depth = to_tune(1L, 3L),
  n.minobsinnode = to_tune(5L, 15L),
  shrinkage = to_tune(0.01, 0.05),

  bag.fraction = 1,  # 关键：每棵树使用全部训练样本
  cv.folds = 0L,     # 内部不再做CV，外层mlr3已经做3折
  n.cores = 1L,
  keep.data = FALSE,
  verbose = FALSE
)

learner_gbm_1 <- ppl(
  "distrcompositor",
  learner = learner_gbm_0,
  estimator = "kaplan",
  form = "ph"
)

learner_gbm_2 <- as_learner(learner_gbm_1)
learner_gbm <- auto_tuner(
  tuner = tnr(
    "grid_search",
    resolution = 2,
    batch_size = 1
  ),
  learner = learner_gbm_2,
  resampling = rsmp("cv", folds = 3),
  measure = msr("surv.cindex"),
  terminator = trm("none")
)

future::plan(future::sequential)

set.seed(seednum)

system.time({
  learner_gbm$train(task_train_gbm)
})

learner_gbm$tuning_result
learner_gbm$archive$data

predtrain_gbm <- learner_gbm$predict(task_train_gbm)
predprobtrain_gbm <-
  predtrain_gbm$distr[1:nrow(traindata2_gbm)]$survival(itps) %>%
  t() %>%
  as.data.frame() %>%
  mutate(model = "gbm",
         dataset = "train",
         time = traindata2_gbm$OS,
         status = traindata2_gbm$status_os)

evaltrain_gbm <- eval4sa(
  predprob = predprobtrain_gbm,
  preddata = traindata2_gbm,
  etime = "OS",
  estatus = "status_os",
  model = "gbm",
  dataset = "train",
  timepoints = itps,
  plotcalimethod = "quantile",
  bw4nne = NULL,
  q4quantile = 1
)

predvalidation2_gbm <- learner_gbm$predict(task_validation2_gbm)
predprobvalidation2_gbm <-
  predvalidation2_gbm$distr[1:nrow(validation_data2_gbm)]$survival(itps) %>%
  t() %>%
  as.data.frame() %>%
  mutate(model = "gbm",
         dataset = "validation2",
         time = validation_data2_gbm$OS,
         status = validation_data2_gbm$status_os)

evalvalidation2_gbm <- eval4sa(
  predprob = predprobvalidation2_gbm,
  preddata = validation_data2_gbm,
  etime = "OS",
  estatus = "status_os",
  model = "gbm",
  dataset = "validation2",
  timepoints = itps,
  plotcalimethod = "quantile",
  bw4nne = NULL,
  q4quantile = 1
)

evaltrain_rsf$auc
evalvalidation2_rsf$auc

save(predtrain_gbm, predprobtrain_gbm, evaltrain_gbm,
     predvalidation2_gbm, predprobvalidation2_gbm, evalvalidation2_gbm,
     file = file.path(model_dir, "gbm.RData"))

cat("========== GBM 模型训练完成 ==========\n\n")

cat("\n========== 开始训练 SVM 模型 ==========\n")

datarecipe_svm <- recipe(OS + status_os ~ ., traindata) %>%
  step_dummy(all_nominal_predictors()) %>%
  step_normalize(all_predictors()) %>%
  prep()

traindata2_svm <- bake(datarecipe_svm, new_data = NULL) %>%
  dplyr::select(OS, status_os, everything())
validation_data2_svm <- bake(datarecipe_svm, new_data = validation_data2_model_data) %>%
  dplyr::select(OS, status_os, everything())

task_train_svm <- as_task_surv(
  traindata2_svm,
  time = "OS",
  event = "status_os",
  type = "right"
)
task_validation2_svm <- as_task_surv(
  validation_data2_svm,
  time = "OS",
  event = "status_os",
  type = "right"
)

learner_svm_0 <- lrn(
  "surv.svm",
  type = "vanbelle1",
  diff.meth = "makediff1",
  kernel = "rbf_kernel",
  opt.meth = "ipop",
  gamma = to_tune(p_dbl(lower = 1e-2, upper = 1))
)

learner_svm_1 <- ppl(
  "distrcompositor",
  learner = learner_svm_0,
  estimator = "kaplan",
  form = "ph"
)

learner_svm_2 <- as_learner(learner_svm_1)

learner_svm <- auto_tuner(
  tuner = tnr("random_search", batch_size = 4),
  learner = learner_svm_2,
  resampling = rsmp("cv", folds = 5),
  measure = msr("surv.cindex"),
  terminator = trm("evals", n_evals = 10)
)

future::plan(future::sequential)

set.seed(seednum)
learner_svm$train(task_train_svm)
learner_svm$tuning_result

predtrain_svm <- learner_svm$predict(task_train_svm)
svm_calibration_data <- data.frame(
  OS = traindata2_svm$OS,
  status_os = traindata2_svm$status_os,
  score = as.numeric(predtrain_svm$crank)
)
if (any(!is.finite(svm_calibration_data$score))) {
  stop("SVM produced non-finite training scores.")
}
# Additional PH calibration fitted only on the development training cohort.
svm_calibration <- survival::coxph(
  survival::Surv(OS, status_os) ~ score,
  data = svm_calibration_data, x = TRUE, y = TRUE, model = TRUE
)
if (any(!is.finite(stats::coef(svm_calibration)))) {
  stop("SVM Cox calibration failed: check score variation and event counts.")
}
svm_baseline_hazard <- survival::basehaz(svm_calibration, centered = TRUE)

predict_svm_probabilities <- function(prediction, timepoints) {
  if (!length(timepoints) || any(!is.finite(timepoints)) || any(timepoints < 0)) {
    stop("Evaluation times must be finite and non-negative.")
  }
  scores <- as.numeric(prediction$crank)
  if (any(!is.finite(scores))) stop("SVM produced non-finite prediction scores.")
  lp <- stats::predict(svm_calibration, newdata = data.frame(score = scores),
                       type = "lp", reference = "sample")
  index <- findInterval(timepoints, svm_baseline_hazard$time)
  hazard <- c(0, svm_baseline_hazard$hazard)[index + 1L]
  probabilities <- exp(-outer(exp(lp), hazard))
  if (any(!is.finite(probabilities))) stop("Non-finite calibrated SVM probabilities.")
  colnames(probabilities) <- as.character(timepoints)
  as.data.frame(probabilities, check.names = FALSE)
}

predprobtrain_svm <-
  predict_svm_probabilities(predtrain_svm, itps) %>%
  mutate(model = "svm",
         dataset = "train",
         time = traindata2_svm$OS,
         status = traindata2_svm$status_os)

evaltrain_svm <- eval4sa(
  predprob = predprobtrain_svm,
  preddata = traindata2_svm,
  etime = "OS",
  estatus = "status_os",
  model = "svm",
  dataset = "train",
  timepoints = itps,
  plotcalimethod = "quantile",
  bw4nne = NULL,
  q4quantile = 1
)

predvalidation2_svm <- learner_svm$predict(task_validation2_svm)
predprobvalidation2_svm <-
  predict_svm_probabilities(predvalidation2_svm, itps) %>%
  mutate(model = "svm",
         dataset = "validation2",
         time = validation_data2_svm$OS,
         status = validation_data2_svm$status_os)

evalvalidation2_svm <- eval4sa(
  predprob = predprobvalidation2_svm,
  preddata = validation_data2_svm,
  etime = "OS",
  estatus = "status_os",
  model = "svm",
  dataset = "validation2",
  timepoints = itps,
  plotcalimethod = "quantile",
  bw4nne = NULL,
  q4quantile = 1
)
evaltrain_svm$auc
evalvalidation2_svm$auc
save(predtrain_svm, predprobtrain_svm, evaltrain_svm,
     predvalidation2_svm, predprobvalidation2_svm, evalvalidation2_svm,
     svm_calibration, svm_baseline_hazard,
     file = file.path(model_dir, "svm.RData"))

cat("========== SVM 模型训练完成 ==========\n\n")

cat("\n========== 开始训练 XGBoost 模型 ==========\n")

datarecipe_xgboost <- recipe(OS + status_os ~ ., traindata) %>%
  step_dummy(all_nominal_predictors()) %>%
  prep()

traindata2_xgboost <- bake(datarecipe_xgboost, new_data = NULL) %>%
  dplyr::select(OS, status_os, everything())
validation_data2_xgboost <- bake(datarecipe_xgboost, new_data = validation_data2_model_data) %>%
  dplyr::select(OS, status_os, everything())

task_train_xgboost <- as_task_surv(
  traindata2_xgboost,
  time = "OS",
  event = "status_os",
  type = "right"
)
task_validation2_xgboost <- as_task_surv(
  validation_data2_xgboost,
  time = "OS",
  event = "status_os",
  type = "right"
)

learner_xgboost_0 <- lrn(
  "surv.xgboost.cox",
  booster = "gbtree",
  nrounds = to_tune(100L, 500L),
  max_depth = to_tune(1L, 5L),
  learning_rate = to_tune(1e-4, 1)
)

learner_xgboost_1 <- ppl(
  "distrcompositor",
  learner = learner_xgboost_0,
  estimator = "kaplan",
  form = "ph"
)

learner_xgboost_2 <- as_learner(learner_xgboost_1)

learner_xgboost <- auto_tuner(
  tuner = tnr("random_search", batch_size = 4),
  learner = learner_xgboost_2,
  resampling = rsmp("cv", folds = 5),
  measure = msr("surv.cindex"),
  terminator = trm("evals", n_evals = 40)
)

future::plan(future::sequential)

set.seed(seednum)
learner_xgboost$train(task_train_xgboost)
learner_xgboost$tuning_result

predtrain_xgboost <- learner_xgboost$predict(task_train_xgboost)
predprobtrain_xgboost <-
  extract_survival_probabilities(predtrain_xgboost, itps, nrow(traindata2_xgboost)) %>%
  mutate(model = "xgboost",
         dataset = "train",
         time = traindata2_xgboost$OS,
         status = traindata2_xgboost$status_os)

evaltrain_xgboost <- eval4sa(
  predprob = predprobtrain_xgboost,
  preddata = traindata2_xgboost,
  etime = "OS",
  estatus = "status_os",
  model = "xgboost",
  dataset = "train",
  timepoints = itps,
  plotcalimethod = "quantile",
  bw4nne = NULL,
  q4quantile = 1
)

predvalidation2_xgboost <- learner_xgboost$predict(task_validation2_xgboost)
predprobvalidation2_xgboost <-
  extract_survival_probabilities(predvalidation2_xgboost, itps, nrow(validation_data2_xgboost)) %>%
  mutate(model = "xgboost",
         dataset = "validation2",
         time = validation_data2_xgboost$OS,
         status = validation_data2_xgboost$status_os)

evalvalidation2_xgboost <- eval4sa(
  predprob = predprobvalidation2_xgboost,
  preddata = validation_data2_xgboost,
  etime = "OS",
  estatus = "status_os",
  model = "xgboost",
  dataset = "validation2",
  timepoints = itps,
  plotcalimethod = "quantile",
  bw4nne = NULL,
  q4quantile = 1
)
evaltrain_xgboost$auc
evalvalidation2_xgboost$auc
save(
  predtrain_xgboost, predprobtrain_xgboost, evaltrain_xgboost,
  predvalidation2_xgboost, predprobvalidation2_xgboost, evalvalidation2_xgboost,
  file = file.path(model_dir, "xgboost.RData"))

validation1_dir <- file.path(output_dir, "validation_set_1")
dir.create(validation1_dir, recursive = TRUE, showWarnings = FALSE)

evaluate_validation1 <- function(
    learner, recipe_object, model_id, raw_validation_data,
    timepoints = itps
) {
  baked_data <- bake(recipe_object, new_data = raw_validation_data) %>%
    dplyr::select(OS, status_os, everything())

  finite_followup <- baked_data$OS[is.finite(baked_data$OS)]
  if (length(finite_followup) == 0L) {
    stop("Validation set 1 has no finite OS follow-up values.")
  }
  max_followup <- max(finite_followup)
  evaluable_timepoints <- sort(unique(
    timepoints[is.finite(timepoints) & timepoints < max_followup]
  ))
  if (length(evaluable_timepoints) == 0L) {
    stop(
      "No requested evaluation time is below the validation set 1 maximum follow-up (",
      round(max_followup, 3), " months)."
    )
  }
  omitted_timepoints <- setdiff(timepoints, evaluable_timepoints)
  if (length(omitted_timepoints) > 0L) {
    message(
      model_id, ": validation set 1 follow-up max = ", round(max_followup, 3),
      " months; evaluating ", paste(evaluable_timepoints, collapse = ", "),
      " months and omitting ", paste(omitted_timepoints, collapse = ", "),
      " months."
    )
  }

  validation_task <- as_task_surv(
    baked_data,
    time = "OS",
    event = "status_os",
    type = "right"
  )
  prediction <- if (identical(model_id, "rpart")) {
    compose_rpart_prediction(validation_task)
  } else {
    learner$predict(validation_task)
  }
  probability <- if (identical(model_id, "svm")) {
    predict_svm_probabilities(prediction, evaluable_timepoints)
  } else {
    if (is.null(prediction$distr)) {
      stop("Model '", model_id, "' returned no survival distribution; predict types: ",
           paste(prediction$predict_types, collapse = ", "))
    }
    extract_survival_probabilities(prediction, evaluable_timepoints, nrow(baked_data))
  }
  probability <- probability %>%
    mutate(
      model = model_id,
      dataset = "validation1",
      time = baked_data$OS,
      status_os = baked_data$status_os
    )
  colnames(probability)[seq_along(evaluable_timepoints)] <- as.character(evaluable_timepoints)

  evaluation <- eval4sa(
    predprob = probability,
    preddata = baked_data,
    etime = "OS",
    estatus = "status_os",
    model = model_id,
    dataset = "validation1",
    timepoints = evaluable_timepoints,
    plotcalimethod = "quantile",
    bw4nne = NULL,
    q4quantile = 3
  )

  list(
    data = baked_data,
    task = validation_task,
    prediction = prediction,
    probability = probability,
    evaluation = evaluation,
    requested_timepoints = timepoints,
    evaluable_timepoints = evaluable_timepoints,
    omitted_timepoints = omitted_timepoints,
    max_followup = max_followup
  )
}

validation1_results <- list(
  coxph = evaluate_validation1(
    learner_coxph, datarecipe_coxph, "coxph", validation_data1
  ),
  lasso = evaluate_validation1(
    learner_lasso, datarecipe_lasso, "lasso", validation_data1
  ),
  rpart = evaluate_validation1(
    learner_rpart, datarecipe_rpart, "rpart", validation_data1
  ),
  rsf = evaluate_validation1(
    learner_rsf, datarecipe_rsf, "rsf", validation_data1
  ),
  gbm = evaluate_validation1(
    learner_gbm, datarecipe_gbm, "gbm", validation_data1
  ),
  svm = evaluate_validation1(
    learner_svm, datarecipe_svm, "svm", validation_data1
  ),
  xgboost = evaluate_validation1(
    learner_xgboost, datarecipe_xgboost, "xgboost", validation_data1
  )
)

validation_data1_coxph <- validation1_results$coxph$data
validation_data1_gbm <- validation1_results$gbm$data
validation_data1_xgboost <- validation1_results$xgboost$data
predvalidation1_coxph <- validation1_results$coxph$prediction
predvalidation1_gbm <- validation1_results$gbm$prediction
predvalidation1_xgboost <- validation1_results$xgboost$prediction
evalvalidation1_coxph <- validation1_results$coxph$evaluation
evalvalidation1_lasso <- validation1_results$lasso$evaluation
evalvalidation1_rpart <- validation1_results$rpart$evaluation
evalvalidation1_rsf <- validation1_results$rsf$evaluation
evalvalidation1_gbm <- validation1_results$gbm$evaluation
evalvalidation1_svm <- validation1_results$svm$evaluation
evalvalidation1_xgboost <- validation1_results$xgboost$evaluation

validation1_auc <- bind_rows(lapply(
  validation1_results,
  function(x) x$evaluation$auc
))
write.csv(
  validation1_auc,
  file.path(validation1_dir, "validation1_AUC.csv"),
  row.names = FALSE
)

for (model_id in names(validation1_results)) {
  model_result <- validation1_results[[model_id]]
  write.csv(
    model_result$probability,
    file.path(validation1_dir, paste0(model_id, "_survival_probability.csv")),
    row.names = FALSE
  )
  write.csv(
    model_result$evaluation$brierscore,
    file.path(validation1_dir, paste0(model_id, "_brier.csv")),
    row.names = FALSE
  )
  ggsave(
    file.path(validation1_dir, paste0(model_id, "_ROC.pdf")),
    model_result$evaluation$rocplot,
    width = 7,
    height = 6
  )
  ggsave(
    file.path(validation1_dir, paste0(model_id, "_calibration.pdf")),
    model_result$evaluation$calibrationplot,
    width = 7,
    height = 6
  )
}

saveRDS(
  validation1_results,
  file.path(validation1_dir, "validation1_all_models.rds")
)

explain_dir <- file.path(output_dir, "explainability")
dir.create(explain_dir, recursive = TRUE, showWarnings = FALSE)

cox_summary <- summary(learner_coxph$model)
cox_explanation <- data.frame(
  variable = rownames(cox_summary$coefficients),
  beta = cox_summary$coefficients[, "coef"],
  HR = cox_summary$conf.int[, "exp(coef)"],
  CI_lower = cox_summary$conf.int[, "lower .95"],
  CI_upper = cox_summary$conf.int[, "upper .95"],
  p_value = cox_summary$coefficients[, "Pr(>|z|)"],
  row.names = NULL,
  check.names = FALSE
) %>%
  arrange(p_value)
write.csv(
  cox_explanation,
  file.path(explain_dir, "coxph_HR_95CI.csv"),
  row.names = FALSE
)

permutation_importance_surv <- function(
    learner, data, time_col = "OS", event_col = "status_os",
    model_name, dataset_name, n_repeats = 20L, seed = 26L
) {
  feature_names <- setdiff(names(data), c(time_col, event_col))
  base_task <- as_task_surv(
    data, time = time_col, event = event_col, type = "right"
  )
  cindex_measure <- msr("surv.cindex")
  baseline_cindex <- unname(
    learner$predict(base_task)$score(cindex_measure)[[1L]]
  )

  set.seed(seed)
  result <- lapply(feature_names, function(feature_name) {
    permuted_cindex <- replicate(n_repeats, {
      permuted_data <- data
      permuted_data[[feature_name]] <- sample(permuted_data[[feature_name]])
      permuted_task <- as_task_surv(
        permuted_data, time = time_col, event = event_col, type = "right"
      )
      unname(learner$predict(permuted_task)$score(cindex_measure)[[1L]])
    })

    data.frame(
      model = model_name,
      dataset = dataset_name,
      variable = feature_name,
      baseline_cindex = baseline_cindex,
      permuted_cindex = mean(permuted_cindex, na.rm = TRUE),
      importance = baseline_cindex - mean(permuted_cindex, na.rm = TRUE),
      importance_sd = stats::sd(
        baseline_cindex - permuted_cindex, na.rm = TRUE
      ),
      stringsAsFactors = FALSE
    )
  })

  bind_rows(result) %>% arrange(desc(importance))
}

permutation_importance_validation1 <- bind_rows(
  permutation_importance_surv(
    learner_coxph, validation_data1_coxph,
    model_name = "CoxPH", dataset_name = "validation1", seed = seednum
  ),
  permutation_importance_surv(
    learner_gbm, validation_data1_gbm,
    model_name = "GBM", dataset_name = "validation1", seed = seednum
  ),
  permutation_importance_surv(
    learner_xgboost, validation_data1_xgboost,
    model_name = "XGBoost", dataset_name = "validation1", seed = seednum
  )
)

write.csv(
  permutation_importance_validation1,
  file.path(explain_dir, "permutation_importance_validation1.csv"),
  row.names = FALSE
)

p_permutation_importance <- permutation_importance_validation1 %>%
  mutate(
    variable = forcats::fct_reorder(variable, importance),
    model = factor(model, levels = c("CoxPH", "GBM", "XGBoost"))
  ) %>%
  ggplot(aes(x = importance, y = variable, color = model)) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey60") +
  geom_point(position = position_dodge(width = 0.6), size = 2.5) +
  geom_errorbar(
    aes(
      xmin = importance - importance_sd,
      xmax = importance + importance_sd
    ),
    position = position_dodge(width = 0.6),
    width = 0
  ) +
  labs(
    x = "Permutation importance (decrease in validation set 1 C-index)",
    y = NULL,
    color = "Model",
    title = "Model-agnostic variable importance"
  ) +
  theme_bw()

ggsave(
  file.path(explain_dir, "permutation_importance_validation1.pdf"),
  p_permutation_importance,
  width = 8,
  height = max(
    4,
    0.35 * length(unique(permutation_importance_validation1$variable)) + 2
  )
)

save(
  cox_explanation,
  permutation_importance_validation1,
  p_permutation_importance,
  file = file.path(explain_dir, "model_explainability.RData")
)

} # End training and analyses preceding SHAP.

if (SHAP_ONLY) {
  required_objects <- c("traindata", "validation_data2_model_data", "seednum",
    paste0("learner_", c("coxph", "lasso", "rpart", "rsf", "gbm", "svm", "xgboost")),
    paste0("datarecipe_", c("coxph", "lasso", "rpart", "rsf", "gbm", "svm", "xgboost")))
  absent <- required_objects[!vapply(required_objects, exists, logical(1), inherits = TRUE)]
  if (length(absent)) stop("SHAP_ONLY requires trained objects: ", paste(absent, collapse = ", "),
                          ". Set SHAP_ONLY <- FALSE to run the full analysis.")
  library(dplyr)
  library(ggplot2)
  explain_dir <- file.path(output_dir, "explainability")
}

shap_packages <- c("kernelshap", "shapviz", "forcats", "ragg", "svglite")
missing_shap_packages <- shap_packages[
  !vapply(shap_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_shap_packages) > 0L) {
  warning(
    "SHAP analysis skipped. Install required R packages and rerun: ",
    paste(missing_shap_packages, collapse = ", "),
    "\ninstall.packages(c('kernelshap', 'shapviz', 'forcats', 'ragg', 'svglite'))"
  )
} else {
  shap_dir <- file.path(explain_dir, "SHAP_original_variables")
  dir.create(shap_dir, recursive = TRUE, showWarnings = FALSE)

  shap_palette <- c(low = "#4B146E", mid = "#C43C63", high = "#FDB515")
  shap_model_palette <- c(CoxPH = "#3B6FB6", LASSO = "#8E63AC", Rpart = "#A6761D",
                          RSF = "#159D89", GBM = "#ED0000", SVM = "#CC79A7",
                          XGBoost = "#D67A2C")

  theme_shap_publication <- function(base_size = 8.5) {
    ggplot2::theme_classic(base_size = base_size, base_family = "Arial") +
      ggplot2::theme(
        axis.line = ggplot2::element_line(linewidth = 0.45, colour = "black"),
        axis.ticks = ggplot2::element_line(linewidth = 0.4, colour = "black"),
        axis.text = ggplot2::element_text(colour = "black"),
        axis.title = ggplot2::element_text(colour = "black"),
        legend.title = ggplot2::element_text(face = "bold"),
        legend.key.height = grid::unit(14, "mm"),
        plot.title = ggplot2::element_text(face = "bold", size = base_size + 1),
        plot.subtitle = ggplot2::element_text(colour = "#555555"),
        plot.margin = ggplot2::margin(7, 9, 7, 7)
      )
  }

  pretty_feature_label <- function(x) {
    x <- gsub("[._]+", " ", x)
    x <- gsub("\\s+", " ", x)
    trimws(x)
  }

  canonical_feature_label <- function(x) {
    sub(" (Yes|No|Low|High)$", "", pretty_feature_label(x))
  }

  save_shap_plot <- function(plot_object, filename, width = 8, height = 6) {
    filename_base <- tools::file_path_sans_ext(filename)
    tryCatch(
      {
        for (extension in c("pdf", "svg", "tiff", "png")) {
          output_file <- paste0(filename_base, ".", extension)
          save_args <- list(
            filename = output_file, plot = plot_object,
            width = width, height = height, units = "in",
            limitsize = FALSE, bg = "white"
          )
          if (extension %in% c("tiff", "png")) save_args$dpi <- 600
          if (extension == "tiff") save_args$compression <- "lzw"
          do.call(ggplot2::ggsave, save_args)
        }
        TRUE
      },
      error = function(e) {
        warning("Could not save SHAP plot ", filename, ": ", conditionMessage(e))
        FALSE
      }
    )
  }

  safe_filename <- function(x) {
    x <- gsub("[^[:alnum:]_.-]+", "_", x)
    substr(x, 1L, 120L)
  }

  mlr3_crank_predict <- function(object, X) {
    # Recreate verified duplicate inputs for already-trained models.
    if (!is.null(object$input_map)) {
      X <- as.data.frame(X)[, unname(object$input_map), drop = FALSE]
      names(X) <- names(object$input_map)
    }
    baked <- recipes::bake(object$recipe, new_data = as.data.frame(X),
                           recipes::all_predictors())
    risk <- as.numeric(object$learner$predict_newdata(as.data.frame(baked))$crank)
    if (length(risk) != nrow(X) || any(!is.finite(risk))) {
      stop("SHAP prediction returned missing or non-finite risk scores.")
    }
    risk
  }

  calculate_and_plot_shap <- function(
    learner, train_data, validation_data2, model_name, recipe_object,
    time_col = "OS", event_col = "status_os",
    background_n = 50L, explain_n = 100L, seed = 26L
  ) {
    model_dir <- file.path(shap_dir, model_name)
    dependence_dir <- file.path(model_dir, "dependence")
    local_dir <- file.path(model_dir, "local")
    dir.create(dependence_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(local_dir, recursive = TRUE, showWarnings = FALSE)

    # Use recipe input names, excluding scores appended by later analyses.
    predictor_names <- unique(recipe_object$var_info$variable[
      recipe_object$var_info$role == "predictor"
    ])
    if (!length(predictor_names) ||
        !all(predictor_names %in% names(train_data)) ||
        !all(predictor_names %in% names(validation_data2))) {
      stop("Original recipe predictors are missing from SHAP input data.")
    }
    input_map <- stats::setNames(predictor_names, predictor_names)
    for (variable in predictor_names) {
      original <- sub("\\.[0-9]+$", "", variable)
      if (original != variable && original %in% predictor_names) {
        same_train <- identical(train_data[[variable]], train_data[[original]])
        same_validation2 <- identical(validation_data2[[variable]], validation_data2[[original]])
        if (!same_train || !same_validation2) {
          stop("Possible duplicate predictor differs from original: ", variable,
               ". Review column provenance before grouping.")
        }
        input_map[[variable]] <- original
      }
    }
    original_names <- unique(unname(input_map))
    train_x <- train_data[, original_names, drop = FALSE]
    test_x <- validation_data2[, original_names, drop = FALSE]
    if (!setequal(names(train_x), names(test_x))) {
      stop("SHAP training and validation set 2 predictor names differ.")
    }
    test_x <- test_x[, names(train_x), drop = FALSE]
    prediction_object <- list(learner = learner, recipe = recipe_object,
                              input_map = input_map)
    message(model_name, " SHAP original variables: ", paste(original_names, collapse = ", "))

    set.seed(seed)
    background_ids <- sample(
      seq_len(nrow(train_x)),
      size = min(background_n, nrow(train_x))
    )
    explain_ids <- if (nrow(test_x) <= explain_n) {
      seq_len(nrow(test_x))
    } else {
      sort(sample(seq_len(nrow(test_x)), explain_n))
    }
    background_x <- train_x[background_ids, , drop = FALSE]
    explain_x <- test_x[explain_ids, , drop = FALSE]

    set.seed(seed)
    shap_kernel <- kernelshap::kernelshap(
      object = prediction_object,
      X = explain_x,
      bg_X = background_x,
      pred_fun = mlr3_crank_predict,
      verbose = FALSE
    )
    shap_object <- shapviz::shapviz(shap_kernel)

    shap_values <- as.data.frame(shap_object$S)
    shap_values <- data.frame(
      explained_row = explain_ids,
      shap_values,
      check.names = FALSE
    )
    write.csv(
      shap_values,
      file.path(model_dir, "SHAP_values.csv"),
      row.names = FALSE
    )

    mean_abs_shap <- data.frame(
      model = model_name,
      variable = colnames(shap_object$S),
      mean_abs_shap = colMeans(abs(shap_object$S), na.rm = TRUE),
      row.names = NULL
    ) %>%
      mutate(
        importance_share = if (sum(mean_abs_shap, na.rm = TRUE) > 0)
          mean_abs_shap / sum(mean_abs_shap, na.rm = TRUE) else 0
      ) %>%
      arrange(desc(mean_abs_shap))
    write.csv(
      mean_abs_shap,
      file.path(model_dir, "SHAP_mean_absolute_importance.csv"),
      row.names = FALSE
    )

    feature_order <- mean_abs_shap$variable
    bar_data <- mean_abs_shap %>%
      mutate(
        variable_label = pretty_feature_label(variable),
        variable_label = factor(variable_label, levels = rev(pretty_feature_label(feature_order)))
      )
    p_bar <- ggplot2::ggplot(
      bar_data,
      ggplot2::aes(x = mean_abs_shap, y = variable_label)
    ) +
      ggplot2::geom_col(width = 0.68, fill = shap_model_palette[[model_name]]) +
      ggplot2::geom_text(
        ggplot2::aes(label = formatC(mean_abs_shap, digits = 3, format = "f")),
        hjust = 1.08, colour = "white", size = 2.8, family = "Arial"
      ) +
      ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.08))) +
      ggplot2::labs(
        x = expression("Mean |SHAP value|"), y = NULL,
        title = paste0(model_name, " — global feature importance"),
        subtitle = "Validation set 2 cohort"
      ) +
      theme_shap_publication()
    save_shap_plot(
      p_bar,
      file.path(model_dir, "SHAP_importance_bar.pdf"),
      width = 6.2,
      height = max(4, 0.27 * ncol(explain_x) + 1.8)
    )

    set.seed(seed)
    beeswarm_display_order <- rev(feature_order)
    beeswarm_rows <- lapply(seq_along(beeswarm_display_order), function(i) {
      feature_name <- beeswarm_display_order[[i]]
      raw_value <- explain_x[[feature_name]]
      numeric_value <- if (is.numeric(raw_value)) {
        as.numeric(raw_value)
      } else {
        as.numeric(factor(raw_value, levels = unique(raw_value)))
      }
      value_range <- range(numeric_value, na.rm = TRUE)
      normalized_value <- if (!all(is.finite(value_range)) || diff(value_range) == 0) {
        rep(0.5, length(numeric_value))
      } else {
        (numeric_value - value_range[1]) / diff(value_range)
      }
      data.frame(
        variable = feature_name,
        variable_label = pretty_feature_label(feature_name),
        shap_value = as.numeric(shap_object$S[, feature_name]),
        feature_value = normalized_value,
        y = i + stats::runif(length(numeric_value), -0.22, 0.22),
        stringsAsFactors = FALSE
      )
    })
    beeswarm_data <- dplyr::bind_rows(beeswarm_rows)
    beeswarm_data$variable_label <- factor(
      beeswarm_data$variable_label,
      levels = pretty_feature_label(beeswarm_display_order)
    )
    beeswarm_data$y <- as.numeric(beeswarm_data$variable_label) +
      (beeswarm_data$y - round(beeswarm_data$y))

    p_beeswarm <- ggplot2::ggplot(
      beeswarm_data,
      ggplot2::aes(x = shap_value, y = y, colour = feature_value)
    ) +
      ggplot2::geom_vline(xintercept = 0, colour = "#A9A9A9", linewidth = 0.55) +
      ggplot2::geom_point(size = 1.55, alpha = 0.82, stroke = 0) +
      ggplot2::scale_y_continuous(
        breaks = seq_along(beeswarm_display_order),
        labels = pretty_feature_label(beeswarm_display_order),
        expand = ggplot2::expansion(add = 0.6)
      ) +
      ggplot2::scale_colour_gradientn(
        colours = unname(shap_palette), limits = c(0, 1),
        breaks = c(0, 1), labels = c("Low", "High"),
        name = "Feature value"
      ) +
      ggplot2::labs(
        x = "SHAP value (impact on predicted risk)", y = NULL,
        title = paste0(model_name, " — SHAP summary"),
        subtitle = "Positive values increase predicted risk; negative values decrease it"
      ) +
      theme_shap_publication() +
      ggplot2::theme(legend.position = "right")
    save_shap_plot(
      p_beeswarm,
      file.path(model_dir, "SHAP_beeswarm.pdf"),
      width = 7.2,
      height = max(4.2, 0.28 * ncol(explain_x) + 2)
    )

    for (feature_name in colnames(explain_x)) {
      feature_file <- safe_filename(feature_name)
      p_dependence <- shapviz::sv_dependence(
        shap_object,
        v = feature_name,
        color_var = NULL
      ) +
        ggplot2::geom_hline(yintercept = 0, colour = "#A9A9A9", linewidth = 0.45) +
        ggplot2::labs(
          title = paste0(model_name, " — ", pretty_feature_label(feature_name)),
          subtitle = "SHAP dependence in the validation set 2 cohort",
          x = pretty_feature_label(feature_name), y = "SHAP value"
        ) +
        theme_shap_publication()
      save_shap_plot(
        p_dependence,
        file.path(dependence_dir, paste0(feature_file, "_dependence.pdf"))
      )

      if (ncol(explain_x) > 1L) {
        p_dependence_interaction <- tryCatch(
          shapviz::sv_dependence(
            shap_object,
            v = feature_name,
            color_var = "auto"
          ) +
            ggplot2::geom_hline(yintercept = 0, colour = "#A9A9A9", linewidth = 0.45) +
            ggplot2::labs(
              title = paste0(model_name, " — ", pretty_feature_label(feature_name)),
              subtitle = "SHAP dependence; colour shows the strongest interaction",
              x = pretty_feature_label(feature_name), y = "SHAP value",
              colour = "Interaction value"
            ) +
            ggplot2::scale_colour_gradientn(colours = unname(shap_palette)) +
            theme_shap_publication(),
          error = function(e) NULL
        )
        if (!is.null(p_dependence_interaction)) {
          save_shap_plot(
            p_dependence_interaction,
            file.path(
              dependence_dir,
              paste0(feature_file, "_dependence_interaction.pdf")
            )
          )
        }
      }
    }

    explained_risk <- mlr3_crank_predict(prediction_object, explain_x)
    local_rows <- c(
      low_risk = which.min(explained_risk),
      median_risk = which.min(abs(explained_risk - stats::median(explained_risk))),
      high_risk = which.max(explained_risk)
    )

    local_index <- data.frame(
      model = model_name,
      explanation = names(local_rows),
      shap_row = as.integer(local_rows),
      validation2_row = explain_ids[as.integer(local_rows)],
      predicted_crank = explained_risk[as.integer(local_rows)],
      row.names = NULL
    )
    write.csv(
      local_index,
      file.path(local_dir, "local_explanation_rows.csv"),
      row.names = FALSE
    )

    for (local_name in names(local_rows)) {
      row_id <- as.integer(local_rows[[local_name]])
      p_waterfall <- shapviz::sv_waterfall(
        shap_object,
        row_id = row_id,
        max_display = ncol(explain_x)
      ) +
        ggplot2::labs(
          title = paste0(model_name, " — ", gsub("_", " ", local_name), " patient"),
          subtitle = "Local SHAP explanation of predicted risk"
        ) +
        theme_shap_publication()
      save_shap_plot(
        p_waterfall,
        file.path(local_dir, paste0(local_name, "_SHAP_waterfall.pdf"))
      )

      p_force <- tryCatch(
        shapviz::sv_force(
          shap_object,
          row_id = row_id,
          max_display = ncol(explain_x)
        ) +
          ggplot2::labs(
            title = paste0(model_name, " — ", gsub("_", " ", local_name), " patient"),
            subtitle = "Yellow increases risk; purple decreases risk"
          ) +
          theme_shap_publication(),
        error = function(e) NULL
      )
      if (!is.null(p_force)) {
        save_shap_plot(
          p_force,
          file.path(local_dir, paste0(local_name, "_SHAP_force.pdf")),
          width = 10,
          height = 4
        )
      }
    }

    saveRDS(shap_kernel, file.path(model_dir, "kernelshap_result.rds"))
    saveRDS(shap_object, file.path(model_dir, "shapviz_object.rds"))
    saveRDS(list(feature_space = "original_predictors", method = "raw_variable_kernelshap_v2", prediction_type = "crank",
                 variables = names(train_x), input_map = input_map, background_rows = background_ids,
                 explained_rows = explain_ids, seed = seed),
            file.path(model_dir, "SHAP_metadata.rds"))

    list(
      kernelshap = shap_kernel,
      shapviz = shap_object,
      importance = mean_abs_shap,
      local_index = local_index
    )
  }

  shap_results <- list(
    CoxPH = calculate_and_plot_shap(
      learner_coxph, traindata, validation_data2_model_data,
      model_name = "CoxPH", recipe_object = datarecipe_coxph, seed = seednum
    ),
    LASSO = calculate_and_plot_shap(
      learner_lasso, traindata, validation_data2_model_data,
      model_name = "LASSO", recipe_object = datarecipe_lasso, seed = seednum
    ),
    Rpart = calculate_and_plot_shap(
      learner_rpart, traindata, validation_data2_model_data,
      model_name = "Rpart", recipe_object = datarecipe_rpart, seed = seednum
    ),
    RSF = calculate_and_plot_shap(
      learner_rsf, traindata, validation_data2_model_data,
      model_name = "RSF", recipe_object = datarecipe_rsf, seed = seednum
    ),
    SVM = calculate_and_plot_shap(
      learner_svm, traindata, validation_data2_model_data,
      model_name = "SVM", recipe_object = datarecipe_svm, seed = seednum
    ),
    GBM = calculate_and_plot_shap(
      learner_gbm, traindata, validation_data2_model_data,
      model_name = "GBM", recipe_object = datarecipe_gbm, seed = seednum
    ),
    XGBoost = calculate_and_plot_shap(
      learner_xgboost, traindata, validation_data2_model_data,
      model_name = "XGBoost", recipe_object = datarecipe_xgboost, seed = seednum
    )
  )

  shap_importance_all_models <- bind_rows(
    lapply(shap_results, `[[`, "importance")
  )
  write.csv(
    shap_importance_all_models,
    file.path(shap_dir, "SHAP_importance_all_models.csv"),
    row.names = FALSE
  )

  p_shap_comparison <- shap_importance_all_models %>%
    mutate(
      model = factor(model, levels = c("CoxPH", "GBM", "XGBoost")),
      variable = canonical_feature_label(variable),
      variable = forcats::fct_reorder(variable, importance_share)
    ) %>%
    ggplot(aes(x = importance_share, y = variable, color = model)) +
    geom_line(
      aes(group = model),
      position = position_dodge(width = 0.55),
      linewidth = 0.45, alpha = 0.45
    ) +
    geom_point(position = position_dodge(width = 0.55), size = 2.3) +
    scale_color_manual(values = shap_model_palette) +
    labs(
      x = "Within-model share of mean absolute SHAP",
      y = NULL,
      color = "Model",
      title = "SHAP importance comparison",
      subtitle = "Importance is normalized within each model"
    ) +
    theme_shap_publication() +
    theme(legend.position = "top")
  save_shap_plot(
    p_shap_comparison,
    file.path(shap_dir, "SHAP_importance_all_models.pdf"),
    height = max(
      4,
      0.35 * length(unique(shap_importance_all_models$variable)) + 2
    )
  )

  save(
    shap_results,
    shap_importance_all_models,
    p_shap_comparison,
    file = file.path(shap_dir, "SHAP_all_results.RData")
  )
}

cat("========== XGBoost 模型训练完成 ==========\n\n")

evaltrain_coxph$auc
evalvalidation2_coxph$auc
evaltrain_gbm$auc
evalvalidation2_gbm$auc
evaltrain_lasso$auc
evalvalidation2_lasso$auc
evaltrain_rpart$auc
evalvalidation2_rpart$auc
evaltrain_rsf$auc
evalvalidation2_rsf$auc
evaltrain_svm$auc
evalvalidation2_svm$auc
evaltrain_xgboost$auc
evalvalidation2_xgboost$auc

all_auc_list <- list(
  evaltrain_coxph$auc,
  evalvalidation2_coxph$auc,
  evaltrain_gbm$auc,
  evalvalidation2_gbm$auc,
  evaltrain_lasso$auc,
  evalvalidation2_lasso$auc,
  evaltrain_rpart$auc,
  evalvalidation2_rpart$auc,
  evaltrain_rsf$auc,
  evalvalidation2_rsf$auc,
  evaltrain_svm$auc,
  evalvalidation2_svm$auc,
  evaltrain_xgboost$auc,
  evalvalidation2_xgboost$auc,
  evalvalidation1_coxph$auc,
  evalvalidation1_gbm$auc,
  evalvalidation1_lasso$auc,
  evalvalidation1_rpart$auc,
  evalvalidation1_rsf$auc,
  evalvalidation1_svm$auc,
  evalvalidation1_xgboost$auc
)
all_auc_df <- bind_rows(all_auc_list)

auc_wide_df <- copy(all_auc_df)
auc_wide_df[, col_name := paste0(times, "month", dataset)]
auc_wide_df <- dcast(
  auc_wide_df[, .(model, col_name, AUC)],
  model ~ col_name,
  value.var = "AUC"
)
auc_columns <- c(
  "model",
  "12monthtrain", "24monthtrain", "36monthtrain", "48monthtrain", "60monthtrain",
  "12monthtest", "24monthtest", "36monthtest","48monthtest", "60monthtest",
  "12monthvalidation2", "24monthvalidation2", "36monthvalidation2","48monthvalidation2", "60monthvalidation2"
)
for (missing_column in setdiff(auc_columns, names(auc_wide_df))) {
  auc_wide_df[, (missing_column) := NA_real_]
}
auc_wide_df <- auc_wide_df[, ..auc_columns]
numeric_auc_columns <- setdiff(names(auc_wide_df), "model")
auc_wide_df[, (numeric_auc_columns) := lapply(
  .SD, round, digits = 3
), .SDcols = numeric_auc_columns]

evaltrain_gbm$auc
evaltrain_gbm$roc
evaltrain_gbm$rocplot
evaltrain_gbm$brierscore
evaltrain_gbm$brierscoretest
evaltrain_gbm$calibration
evaltrain_gbm$calibrationplot
evalvalidation2_gbm$calibrationplot
evaltrain_coxph$calibrationplot
evalvalidation2_coxph$calibrationplot

traindata$risk_coxph <- predtrain_coxph$crank
traindata$risk_gbm <- predtrain_gbm$crank
traindata$risk_lasso <- predtrain_lasso$crank
traindata$risk_rpart <- predtrain_rpart$crank
traindata$risk_rsf <- predtrain_rsf$crank
traindata$risk_svm <- predtrain_svm$crank
traindata$risk_xgboost <- predtrain_xgboost$crank

validation_data2_model_data$risk_coxph <- predvalidation2_coxph$crank
validation_data2_model_data$risk_gbm <- predvalidation2_gbm$crank
validation_data2_model_data$risk_lasso <- predvalidation2_lasso$crank
validation_data2_model_data$risk_rpart <- predvalidation2_rpart$crank
validation_data2_model_data$risk_rsf <- predvalidation2_rsf$crank
validation_data2_model_data$risk_svm <- predvalidation2_svm$crank
validation_data2_model_data$risk_xgboost <- predvalidation2_xgboost$crank

validation_data1$risk_coxph <- validation1_results$coxph$prediction$crank
validation_data1$risk_gbm <- validation1_results$gbm$prediction$crank
validation_data1$risk_lasso <- validation1_results$lasso$prediction$crank
validation_data1$risk_rpart <- validation1_results$rpart$prediction$crank
validation_data1$risk_rsf <- validation1_results$rsf$prediction$crank
validation_data1$risk_svm <- validation1_results$svm$prediction$crank
validation_data1$risk_xgboost <- validation1_results$xgboost$prediction$crank

time_points <- c(6,12,18, 24,30, 36,42, 48, 54,60)
time_points <- c(12,18, 24,30, 36,42, 48, 54,60)
lancet_colors <- c(
  "#00468B", "#ED0000", "#42B540", "#0099B4", "#925E9F",
  "#FDAF91", "#AD002A"
)

plot_roc_table <- function(data, data_name, requested_time_points = time_points) {
  max_followup <- max(data$OS[is.finite(data$OS)], na.rm = TRUE)
  plot_time_points <- requested_time_points[requested_time_points < max_followup]
  if (length(plot_time_points) == 0L) {
    stop("No ROC evaluation time is below maximal follow-up for ", data_name, ".")
  }
  omitted_time_points <- setdiff(requested_time_points, plot_time_points)
  if (length(omitted_time_points) > 0L) {
    message(
      data_name, ": ROC omitted non-estimable time points: ",
      paste(omitted_time_points, collapse = ", "), " months."
    )
  }
  model_list <- list(
    CoxPH = data$risk_coxph,
    GBM = data$risk_gbm,
    LASSO = data$risk_lasso,
    RPart = data$risk_rpart,
    RSF = data$risk_rsf,
    SVM = data$risk_svm,
    XGBoost = data$risk_xgboost
  )

  auc_ci_list <- lapply(names(model_list), function(model_name) {
    roc_obj <- timeROC(
      T = data$OS, delta = data$status_os,
      marker = model_list[[model_name]], cause = 1,
      times = plot_time_points, iid = TRUE
    )

    ci_obj <- confint(roc_obj)
    auc_vals <- roc_obj$AUC
    ci_vals <- if (!is.null(ci_obj$CI_AUC)) ci_obj$CI_AUC else ci_obj$CB_AUC

    df_temp <- data.frame(
      Time = plot_time_points,
      Model = model_name,
      AUC_raw = auc_vals,
      AUC_lower = round(ci_vals[, "2.5%"] / 100, 3),
      AUC_upper = round(ci_vals[, "97.5%"] / 100, 3)
    ) %>%
      mutate(
        AUC = round(AUC_raw, 3),
        `AUC (95%CI)` = paste0(AUC, " (", AUC_lower, "–", AUC_upper, ")"),
        Time_label = paste0(Time, " months")
      )
    return(df_temp)
  })

  auc_ci_combined <- do.call(rbind, auc_ci_list)
  auc_ci_results_plot <- auc_ci_combined %>% dplyr::select(Model, Time, AUC)
  auc_ci_results_table <- auc_ci_combined %>%
    dplyr::select(Model, Time_label, `AUC (95%CI)`) %>%
    pivot_wider(names_from = Time_label, values_from = `AUC (95%CI)`)

  roc_plot <- ggplot(auc_ci_results_plot, aes(x = Time, y = AUC, color = Model)) +
    geom_line(linewidth = 1.2, alpha = 0.8) +
    geom_point(size = 2.5, alpha = 0.9) +
    scale_color_manual(values = lancet_colors) +
    scale_y_continuous(limits = c(0.5, 1.0), breaks = seq(0.5, 1.0, 0.1), expand = c(0, 0)) +
    scale_x_continuous(
      limits = range(plot_time_points) + c(-1, 1),
      breaks = plot_time_points,
      expand = c(0, 0)
    ) +
    theme_bw() +
    theme(
      text = element_text(family = "sans", size = 10),
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.direction = "horizontal",
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(color = "black"),
      panel.background = element_rect(fill = "white"),
      plot.margin = unit(c(0.5, 0.5, 0, 0.5), "cm")
    ) +
    labs(
      title = paste0("Time-dependent AUC of Multiple Prognostic Models (", data_name, ")"),
      x = "Time (months)",
      y = "Area Under the ROC Curve (AUC)"
    )

  lancet_table_theme <- ttheme_minimal(
    core = list(fg_params = list(fontsize = 9, fontfamily = "sans")),
    colhead = list(fg_params = list(fontsize = 10, fontface = "bold"), bg_params = list(fill = "#f0f0f0")),
    rowhead = list(fg_params = list(fontsize = 10, fontface = "bold"))
  )
  auc_table <- tableGrob(
    auc_ci_results_table,
    theme = lancet_table_theme,
    rows = NULL,
    cols = c("Model", paste0(plot_time_points, " months"))
  )
  auc_table$widths <- unit(rep(1, ncol(auc_table)), "null")

  final_plot <- grid.arrange(roc_plot, auc_table, ncol = 1, heights = c(4, 2))
  return(final_plot)
}

p_train <- plot_roc_table(traindata, "Training Set")
print(p_train)

p_validation2 <- plot_roc_table(validation_data2_model_data, "Validation Set 2")
print(p_validation2)

p_validation1 <- plot_roc_table(validation_data1, "Validation Set 1")
print(p_validation1)

ggsave(file.path(figure_dir, "TRAIN_ROC.pdf"), p_train, width=18, height=10, dpi=300)
ggsave(file.path(figure_dir, "Validation_Set_2_ROC.pdf"), p_validation2, width=18, height=10, dpi=300)
ggsave(
  file.path(figure_dir, "Validation_Set_1_ROC.pdf"),
  p_validation1,
  width = 18,
  height = 10,
  dpi = 300
)

lancet_colors <- c("#00468B", "#ED0000", "#42B540", "#0099B4", "#925E9F", "#FDAF91", "#AD002A")
model_names <- c("CoxPH", "GBM", "LASSO", "RPart", "RSF", "SVM", "XGBoost")
risk_cols <- c("risk_coxph", "risk_gbm", "risk_lasso", "risk_rpart", "risk_rsf", "risk_svm", "risk_xgboost")


model_names <- c("CoxPH", "GBM", "LASSO", "RPart", "RSF", "SVM", "XGBoost")
risk_vars   <- c("risk_coxph","risk_gbm","risk_lasso","risk_rpart","risk_rsf","risk_svm","risk_xgboost")
colors      <- c("#00468B","#ED0000","#42B540","#0099B4","#925E9F","#FDAF91","#AD002A","black","grey")

dca_train_12 <- dca(
  coxph(Surv(OS, status_os) ~ risk_coxph,    data=traindata),
  coxph(Surv(OS, status_os) ~ risk_gbm,      data=traindata),
  coxph(Surv(OS, status_os) ~ risk_lasso,    data=traindata),
  coxph(Surv(OS, status_os) ~ risk_rpart,    data=traindata),
  coxph(Surv(OS, status_os) ~ risk_rsf,      data=traindata),
  coxph(Surv(OS, status_os) ~ risk_svm,      data=traindata),
  coxph(Surv(OS, status_os) ~ risk_xgboost,  data=traindata),
  model.names = model_names, times=12
)
p_train_12 <- ggplot(dca_train_12) +
  scale_color_manual(values=colors) + theme_bw() +
  ggtitle("Training Set - 12 Months")

dca_train_36 <- dca(
  coxph(Surv(OS, status_os) ~ risk_coxph,    data=traindata),
  coxph(Surv(OS, status_os) ~ risk_gbm,      data=traindata),
  coxph(Surv(OS, status_os) ~ risk_lasso,    data=traindata),
  coxph(Surv(OS, status_os) ~ risk_rpart,    data=traindata),
  coxph(Surv(OS, status_os) ~ risk_rsf,      data=traindata),
  coxph(Surv(OS, status_os) ~ risk_svm,      data=traindata),
  coxph(Surv(OS, status_os) ~ risk_xgboost,  data=traindata),
  model.names = model_names, times=36
)
p_train_36 <- ggplot(dca_train_36) +
  scale_color_manual(values=colors) + theme_bw() +
  ggtitle("Training Set - 36 Months")

dca_train_60 <- dca(
  coxph(Surv(OS, status_os) ~ risk_coxph,    data=traindata),
  coxph(Surv(OS, status_os) ~ risk_gbm,      data=traindata),
  coxph(Surv(OS, status_os) ~ risk_lasso,    data=traindata),
  coxph(Surv(OS, status_os) ~ risk_rpart,    data=traindata),
  coxph(Surv(OS, status_os) ~ risk_rsf,      data=traindata),
  coxph(Surv(OS, status_os) ~ risk_svm,      data=traindata),
  coxph(Surv(OS, status_os) ~ risk_xgboost,  data=traindata),
  model.names = model_names, times=60
)
p_train_60 <- ggplot(dca_train_60) +
  scale_color_manual(values=colors) + theme_bw() +
  ggtitle("Training Set - 60 Months")

dca_test_12 <- dca(
  coxph(Surv(OS, status_os) ~ risk_coxph,    data=validation_data2_model_data),
  coxph(Surv(OS, status_os) ~ risk_gbm,      data=validation_data2_model_data),
  coxph(Surv(OS, status_os) ~ risk_lasso,    data=validation_data2_model_data),
  coxph(Surv(OS, status_os) ~ risk_rpart,    data=validation_data2_model_data),
  coxph(Surv(OS, status_os) ~ risk_rsf,      data=validation_data2_model_data),
  coxph(Surv(OS, status_os) ~ risk_svm,      data=validation_data2_model_data),
  coxph(Surv(OS, status_os) ~ risk_xgboost,  data=validation_data2_model_data),
  model.names = model_names, times=12
)
p_test_12 <- ggplot(dca_test_12) +
  scale_color_manual(values=colors) + theme_bw() +
  ggtitle("Validation Set 2 - 12 Months")

dca_test_36 <- dca(
  coxph(Surv(OS, status_os) ~ risk_coxph,    data=validation_data2_model_data),
  coxph(Surv(OS, status_os) ~ risk_gbm,      data=validation_data2_model_data),
  coxph(Surv(OS, status_os) ~ risk_lasso,    data=validation_data2_model_data),
  coxph(Surv(OS, status_os) ~ risk_rpart,    data=validation_data2_model_data),
  coxph(Surv(OS, status_os) ~ risk_rsf,      data=validation_data2_model_data),
  coxph(Surv(OS, status_os) ~ risk_svm,      data=validation_data2_model_data),
  coxph(Surv(OS, status_os) ~ risk_xgboost,  data=validation_data2_model_data),
  model.names = model_names, times=36
)
p_test_36 <- ggplot(dca_test_36) +
  scale_color_manual(values=colors) + theme_bw() +
  ggtitle("Validation Set 2 - 36 Months")

dca_test_60 <- dca(
  coxph(Surv(OS, status_os) ~ risk_coxph,    data=validation_data2_model_data),
  coxph(Surv(OS, status_os) ~ risk_gbm,      data=validation_data2_model_data),
  coxph(Surv(OS, status_os) ~ risk_lasso,    data=validation_data2_model_data),
  coxph(Surv(OS, status_os) ~ risk_rpart,    data=validation_data2_model_data),
  coxph(Surv(OS, status_os) ~ risk_rsf,      data=validation_data2_model_data),
  coxph(Surv(OS, status_os) ~ risk_svm,      data=validation_data2_model_data),
  coxph(Surv(OS, status_os) ~ risk_xgboost,  data=validation_data2_model_data),
  model.names = model_names, times=60
)
p_test_60 <- ggplot(dca_test_60) +
  scale_color_manual(values=colors) + theme_bw() +
  ggtitle("Validation Set 2 - 60 Months")

dca_val_12 <- dca(
  coxph(Surv(OS, status_os) ~ risk_coxph,    data=validation_data1),
  coxph(Surv(OS, status_os) ~ risk_gbm,      data=validation_data1),
  coxph(Surv(OS, status_os) ~ risk_lasso,    data=validation_data1),
  coxph(Surv(OS, status_os) ~ risk_rpart,    data=validation_data1),
  coxph(Surv(OS, status_os) ~ risk_rsf,      data=validation_data1),
  coxph(Surv(OS, status_os) ~ risk_svm,      data=validation_data1),
  coxph(Surv(OS, status_os) ~ risk_xgboost,  data=validation_data1),
  model.names = model_names, times=12
)
p_val_12 <- ggplot(dca_val_12) +
  scale_color_manual(values=colors) + theme_bw() +
  ggtitle("Validation Set 1 - 12 Months")

dca_val_36 <- dca(
  coxph(Surv(OS, status_os) ~ risk_coxph,    data=validation_data1),
  coxph(Surv(OS, status_os) ~ risk_gbm,      data=validation_data1),
  coxph(Surv(OS, status_os) ~ risk_lasso,    data=validation_data1),
  coxph(Surv(OS, status_os) ~ risk_rpart,    data=validation_data1),
  coxph(Surv(OS, status_os) ~ risk_rsf,      data=validation_data1),
  coxph(Surv(OS, status_os) ~ risk_svm,      data=validation_data1),
  coxph(Surv(OS, status_os) ~ risk_xgboost,  data=validation_data1),
  model.names = model_names, times=36
)
p_val_36 <- ggplot(dca_val_36) +
  scale_color_manual(values=colors) + theme_bw() +
  ggtitle("Validation Set 1 - 36 Months")

dca_val_60 <- dca(
  coxph(Surv(OS, status_os) ~ risk_coxph,    data=validation_data1),
  coxph(Surv(OS, status_os) ~ risk_gbm,      data=validation_data1),
  coxph(Surv(OS, status_os) ~ risk_lasso,    data=validation_data1),
  coxph(Surv(OS, status_os) ~ risk_rpart,    data=validation_data1),
  coxph(Surv(OS, status_os) ~ risk_rsf,      data=validation_data1),
  coxph(Surv(OS, status_os) ~ risk_svm,      data=validation_data1),
  coxph(Surv(OS, status_os) ~ risk_xgboost,  data=validation_data1),
  model.names = model_names, times=60
)
p_val_60 <- ggplot(dca_val_60) +
  scale_color_manual(values=colors) + theme_bw() +
  ggtitle("Validation Set 1 - 60 Months")

p_train_12
p_train_36
p_train_60

p_test_12
p_test_36
p_test_60

allplots <- list(
  p_train_12, p_train_36, p_train_60,
  p_val_12,  p_val_36,  p_val_60,
  p_test_12,  p_test_36,  p_test_60
)

ok <- sapply(allplots, inherits, "ggplot")
if(all(ok)){
  combined <- grid.arrange(
    p_train_12, p_train_36, p_train_60,
    p_test_12,  p_test_36,  p_test_60,
    p_val_12,  p_val_36,  p_val_60,
    ncol=3, nrow=3,
    top="Decision Curve Analysis (7 models)"
  )
  ggsave(file.path(figure_dir, "DCA_6plots.pdf"), combined, width=18, height=10, dpi=300)
} else {
  cat("失败的图：", which(!ok), "\n")
}

Charcutpoint <- surv_cutpoint(traindata, time = "OS", event = "status_os",
                              variables = "risk_xgboost")
plot(Charcutpoint)
nomocutoff <- as.numeric(Charcutpoint["risk_xgboost"][[1]]["estimate"][[1]][[1]])

traindata$xgboost_group <- ifelse(traindata$risk_xgboost>=nomocutoff,"High","Low")
validation_data2_model_data$xgboost_group <- ifelse(validation_data2_model_data$risk_xgboost>=nomocutoff,"High","Low")
validation_data1$xgboost_group <- ifelse(validation_data1$risk_xgboost>=nomocutoff,"High","Low")

fit <- survfit(Surv(traindata$OS,traindata$status_os)~xgboost_group,data = traindata)
fit %>% summary(times = c(12,24,36,48,60))
ggsurvplot(fit,linetype = c(1,1),
           risk.table=T,
           data = traindata,
           title="Training Cohort",
           ylab="Mace (percentage)",xlab = "Time (Months)",
           pval = T,
           xlim = c(0,60),
           break.time.by = 12,
           palette = c("#3c5769","#b40000"))
fit <- survfit(Surv(validation_data2_model_data$OS,validation_data2_model_data$status_os)~xgboost_group,data = validation_data2_model_data)
fit %>% summary(times = c(12,24,36,48,60))
ggsurvplot(fit,linetype = c(1,1),
           risk.table=T,
           data = validation_data2_model_data,
           title="Validation Cohort 2",
           ylab="Mace (percentage)",xlab = "Time (Months)",
           pval = T,
           xlim = c(0,60),
           break.time.by = 12,
           palette = c("#3c5769","#b40000"))
fit <- survfit(Surv(validation_data1$OS,validation_data1$status_os)~xgboost_group,data = validation_data1)
fit %>% summary(times = c(12,24,36,48,60))
ggsurvplot(fit,linetype = c(1,1),
           risk.table=T,
           data = validation_data1,
           title="Validation Cohort 1",
           ylab="Mace (percentage)",xlab = "Time (Months)",
           pval = T,
           xlim = c(0,60),
           break.time.by = 12,
           palette = c("#3c5769","#b40000"))



auc.folds <- c("T12","T24","T36","T48","T60","V12","V24","V36","V48","V60")
auc.folds.data <- data.frame(auc.folds)
auc.folds.data$folds1 <- NA
auc.folds.data$folds2 <- NA
auc.folds.data$folds3 <- NA
auc.folds.data$folds4 <- NA
auc.folds.data$folds5 <- NA
auc.folds.data <- data.frame(t(auc.folds.data))
colnames(auc.folds.data) <- auc.folds.data[1,]
auc.folds.data <- auc.folds.data[-1,]
set.seed(40)
K_folds=5
folds.data <- traindata
folds.data$random_num <- sample(1:nrow(folds.data),nrow(folds.data),replace=F)
folds.data$group <- folds.data$random_num %% K_folds+1
model_roc <- data.frame()
for (i in min(folds.data$group):max(folds.data$group)) {
  folds.data.validation <- subset(folds.data,folds.data$group==i)
  folds.data.train <- subset(folds.data,folds.data$group!=i)
  y.train.roc <- timeROC(T = folds.data.train$OS,
                         delta = folds.data.train$status_os,
                         marker = folds.data.train$risk_xgboost,
                         cause = 1,
                         times = c(12,24,36,48,60),
                         iid = T)
  y.validation.roc <- timeROC(T = folds.data.validation$OS,
                              delta = folds.data.validation$status_os,
                              marker = folds.data.validation$risk_xgboost,
                              cause = 1,
                              times = c(12,24,36,48,60),
                              iid = T)
  auc.folds.data[i,1] <- round(y.train.roc$AUC[[1]],3)
  auc.folds.data[i,2] <- round(y.train.roc$AUC[[2]],3)
  auc.folds.data[i,3] <- round(y.train.roc$AUC[[3]],3)
  auc.folds.data[i,4] <- round(y.train.roc$AUC[[4]],3)
  auc.folds.data[i,5] <- round(y.train.roc$AUC[[5]],3)
  auc.folds.data[i,6] <- round(y.validation.roc$AUC[[1]],3)
  auc.folds.data[i,7] <- round(y.validation.roc$AUC[[2]],3)
  auc.folds.data[i,8] <- round(y.validation.roc$AUC[[3]],3)
  auc.folds.data[i,9] <- round(y.validation.roc$AUC[[4]],3)
  auc.folds.data[i,10] <- round(y.validation.roc$AUC[[5]],3)
}
model_roc <- data.frame("T12"=as.numeric(auc.folds.data[,1]),
                        "T24"=as.numeric(auc.folds.data[,2]),
                        "T36"=as.numeric(auc.folds.data[,3]),
                        "T48"=as.numeric(auc.folds.data[,4]),
                        "T60"=as.numeric(auc.folds.data[,5]),
                        "V12"=as.numeric(auc.folds.data[,6]),
                        "V24"=as.numeric(auc.folds.data[,7]),
                        "V36"=as.numeric(auc.folds.data[,8]),
                        "V48"=as.numeric(auc.folds.data[,9]),
                        "V60"=as.numeric(auc.folds.data[,10]))
summary(model_roc)
train.model.roc <- model_roc[,1:5]
names(train.model.roc) <- c("12Months","24Months","36Months","48Months","60Months")
train.model.roc$group <- "Train"
validation_model_roc <- model_roc[,6:10]
names(validation_model_roc) <- c("12Months","24Months","36Months","48Months","60Months")
validation_model_roc$group <- "Validation"
model_roc_plot <- rbind.data.frame(train.model.roc,validation_model_roc)
model_roc_plot <- reshape2::melt(model_roc_plot)
p1<-ggplot(model_roc_plot,aes(group,value,color=group,fill=group))+
  geom_bar(stat="summary",fun=mean,position="dodge")+
  stat_summary(geom = "errorbar",fun.data = 'mean_sd', width = 0.3)+
  labs(x="5 Cross Validation AUC",y=NULL)+
  theme_prism(palette = "candy_bright",
              base_fontface = "bold",
              base_family = "Arial",
              base_size = 16,
              base_line_size = 0.8,
              axis_text_angle = 45)+
  scale_fill_prism(palette = "candy_bright")
p2<-p1+facet_grid(~variable,scales = 'free')
p2
p3 <- p2+geom_point(data=model_roc_plot,aes(group,value),size=2,pch=20,color="black")
p3

dd<-datadist(traindata)
options(datadist = 'dd')
trainnomo<-cph(Surv(OS,status_os)~risk_xgboost,x=T,y=T,data = traindata,surv = T,time.inc = 12)
cal_t_24 <- calibrate(trainnomo, cmethod='KM', method='boot', u=12, m=60, B=1000)
plot(cal_t_24,lwd=2,lty=1,errbar.col=NA,
     xlim=c(0,1),ylim=c(0,1),
     xlab="预测概率",
     ylab="实际发生概率",
     col="#3B83C0")
trainnomo<-cph(Surv(OS,status_os)~risk_xgboost,x=T,y=T,data = traindata,surv = T,time.inc = 24)
cal_t_12 <- calibrate(trainnomo, cmethod='KM', method='boot', u=24, m=60, B=1000)
plot(cal_t_12,add=T,lwd=2,lty=1,xlim=c(0,1),ylim=c(0,1),errbar.col=NA,col="#EE362B")
trainnomo<-cph(Surv(OS,status_os)~risk_xgboost,x=T,y=T,data = traindata,surv = T,time.inc = 36)
cal_t_12 <- calibrate(trainnomo, cmethod='KM', method='boot', u=36, m=60, B=1000)
plot(cal_t_12,add=T,lwd=2,lty=1,xlim=c(0,1),ylim=c(0,1),errbar.col=NA,col="#F8CB75")
legend(0.6,0.4,c(paste("列线图（1年）"),
                 paste("列线图（2年）"),
                 paste("列线图（3年）")),
       x.intersp = 1,
       y.intersp = 1,
       lty = 5,lwd = 3,col = c("#3B83C0","#EE362B","#F8CB75"),bty = "n",
       seg.len = 1,cex = 0.9)
legend("topleft",c(paste("训练集")),x.intersp = 1,
       y.intersp = 0.8,lty = 1,lwd = 2,col = "black",
       bty = "n",seg.len = 1,cex = 0.8)

dd<-datadist(validation_data2_model_data)
options(datadist = 'dd')
trainnomo<-cph(Surv(OS,status_os)~risk_xgboost,x=T,y=T,data = validation_data2_model_data,surv = T,time.inc = 12)
cal_t_24 <- calibrate(trainnomo, cmethod='KM', method='boot', u=12, m=27, B=1000)
plot(cal_t_24,lwd=2,lty=1,errbar.col=NA,
     xlim=c(0,1),ylim=c(0,1),
     xlab="预测概率",
     ylab="实际发生概率",
     col="#3B83C0")
trainnomo<-cph(Surv(OS,status_os)~risk_xgboost,x=T,y=T,data = validation_data2_model_data,surv = T,time.inc = 24)
cal_t_12 <- calibrate(trainnomo, cmethod='KM', method='boot', u=24, m=27, B=1000)
plot(cal_t_12,add=T,lwd=2,lty=1,errbar.col=NA,col="#EE362B")
trainnomo<-cph(Surv(OS,status_os)~risk_xgboost,x=T,y=T,data = validation_data2_model_data,surv = T,time.inc = 36)
cal_t_12 <- calibrate(trainnomo, cmethod='KM', method='boot', u=36, m=27, B=1000)
plot(cal_t_12,add=T,lwd=2,lty=1,errbar.col=NA,col="#F8CB75")
legend(0.6,0.4,c(paste("列线图（1年）"),
                 paste("列线图（2年）"),
                 paste("列线图（3年）")),
       x.intersp = 1,
       y.intersp = 1,
       lty = 5,lwd = 3,col = c("#3B83C0","#EE362B","#F8CB75"),bty = "n",
       seg.len = 1,cex = 0.9)
legend("topleft",c(paste("验证集")),x.intersp = 1,
       y.intersp = 0.8,lty = 1,lwd = 2,col = "black",
       bty = "n",seg.len = 1,cex = 0.8)

dd<-datadist(validation_data1)
options(datadist = 'dd')
trainnomo<-cph(Surv(OS,status_os)~risk_xgboost,x=T,y=T,data = validation_data1,surv = T,time.inc = 12)
cal_t_24 <- calibrate(trainnomo, cmethod='KM', method='boot', u=12, m=27, B=1000)
plot(cal_t_24,lwd=2,lty=1,errbar.col=NA,
     xlim=c(0,1),ylim=c(0,1),
     xlab="预测概率",
     ylab="实际发生概率",
     col="#3B83C0")
trainnomo<-cph(Surv(OS,status_os)~risk_xgboost,x=T,y=T,data = validation_data1,surv = T,time.inc = 24)
cal_t_12 <- calibrate(trainnomo, cmethod='KM', method='boot', u=24, m=27, B=1000)
plot(cal_t_12,add=T,lwd=2,lty=1,errbar.col=NA,col="#EE362B")
trainnomo<-cph(Surv(OS,status_os)~risk_xgboost,x=T,y=T,data = validation_data1,surv = T,time.inc = 36)
cal_t_12 <- calibrate(trainnomo, cmethod='KM', method='boot', u=36, m=27, B=1000)
plot(cal_t_12,add=T,lwd=2,lty=1,errbar.col=NA,col="#F8CB75")

legend(0.6,0.4,c(paste("列线图（1年）"),
                 paste("列线图（2年）"),
                 paste("列线图（3年）")),
       x.intersp = 1,
       y.intersp = 1,
       lty = 5,lwd = 3,col = c("#3B83C0","#EE362B","#F8CB75"),bty = "n",
       seg.len = 1,cex = 0.9)
legend("topleft",c(paste("Validation")),x.intersp = 1,
       y.intersp = 0.8,lty = 1,lwd = 2,col = "black",
       bty = "n",seg.len = 1,cex = 0.8)
# Regenerate publication-ready SHAP figures without retraining models.
required_packages <- c(
  "ggplot2", "dplyr", "forcats", "shapviz", "ragg", "svglite"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Install required packages: ", paste(missing_packages, collapse = ", "))
}
invisible(lapply(required_packages, library, character.only = TRUE))

input_root <- file.path(output_dir, "explainability", "SHAP_original_variables")
output_root <- file.path(output_dir, "explainability", "SHAP_publication_original_variables")
if (!dir.exists(input_root)) {
  stop("SHAP input directory not found: ", input_root)
}
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

models <- c("CoxPH", "LASSO", "Rpart", "RSF", "GBM", "SVM", "XGBoost")
model_colours <- c(CoxPH = "#3B6FB6", LASSO = "#8E63AC", Rpart = "#A6761D",
                   RSF = "#159D89", GBM = "#ED0000", SVM = "#CC79A7", XGBoost = "#D67A2C")
value_colours <- c("#4B146E", "#C43C63", "#FDB515")



pretty_label <- function(x) {
  trimws(gsub("\\s+", " ", gsub("[._]+", " ", x)))
}

canonical_label <- function(x) {
  sub(" (Yes|No|Low|High)$", "", pretty_label(x))
}

theme_shap <- function(base_size = 8.5) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.45, colour = "black"),
      axis.ticks = element_line(linewidth = 0.4, colour = "black"),
      axis.text = element_text(colour = "black"),
      axis.title = element_text(colour = "black"),
      legend.title = element_text(face = "bold"),
      legend.key.height = grid::unit(14, "mm"),
      plot.title = element_text(face = "bold", size = base_size + 1),
      plot.subtitle = element_text(colour = "#555555"),
      plot.margin = margin(7, 9, 7, 7)
    )
}

save_publication <- function(p, base, width, height) {
  dir.create(dirname(base), recursive = TRUE, showWarnings = FALSE)
  ggsave(paste0(base, ".pdf"), p, width = width, height = height,
         units = "in", device = grDevices::cairo_pdf, bg = "white")
  ggsave(paste0(base, ".svg"), p, width = width, height = height,
         units = "in", device = svglite::svglite, bg = "white")
  ragg::agg_tiff(paste0(base, ".tiff"), width = width, height = height,
                 units = "in", res = 600, compression = "lzw")
  print(p)
  dev.off()
  ragg::agg_png(paste0(base, ".png"), width = width, height = height,
                units = "in", res = 300)
  print(p)
  dev.off()
}

importance_list <- list()

for (model_name in models) {
  model_input <- file.path(input_root, model_name)
  model_output <- file.path(output_root, model_name)
  required_files <- file.path(
    model_input,
    c("shapviz_object.rds", "SHAP_mean_absolute_importance.csv")
  )
  if (any(!file.exists(required_files))) {
    stop(
      model_name, " SHAP files are missing: ",
      paste(required_files[!file.exists(required_files)], collapse = ", ")
    )
  }
  dir.create(model_output, recursive = TRUE, showWarnings = FALSE)

  sv <- readRDS(file.path(model_input, "shapviz_object.rds"))
  metadata <- readRDS(file.path(model_input, "SHAP_metadata.rds"))
  if (!identical(metadata$method, "raw_variable_kernelshap_v2") ||
      !identical(colnames(sv$S), metadata$variables)) {
    stop("Recompute SHAP with the updated original-variable code: ", model_name)
  }
  importance <- data.frame(model = model_name, variable = colnames(sv$S),
                           mean_abs_shap = colMeans(abs(sv$S)))
  total_importance <- sum(importance$mean_abs_shap)
  importance$importance_share <- if (total_importance > 0)
    importance$mean_abs_shap / total_importance else 0
  importance <- dplyr::arrange(importance, desc(mean_abs_shap))
  importance_list[[model_name]] <- importance
  order_top <- importance$variable
  order_bottom <- rev(order_top)

  bar_data <- importance %>%
    mutate(label = factor(pretty_label(variable), levels = pretty_label(order_bottom)))
  p_bar <- ggplot(bar_data, aes(mean_abs_shap, label)) +
    geom_col(width = 0.68, fill = model_colours[[model_name]]) +
    geom_text(aes(label = formatC(mean_abs_shap, digits = 3, format = "f")),
              hjust = 1.08, colour = "white", size = 2.8, family = "Arial") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.08))) +
    labs(x = "Mean |SHAP value|", y = NULL,
         title = paste0(model_name, " — global feature importance"),
         subtitle = "Validation set 2 cohort") +
    theme_shap()
  save_publication(p_bar, file.path(model_output, "SHAP_importance"),
                   6.2, max(4, 0.27 * nrow(importance) + 1.8))

  set.seed(26)
  bee <- bind_rows(lapply(seq_along(order_bottom), function(i) {
    feature <- order_bottom[[i]]
    raw <- sv$X[[feature]]
    numeric_raw <- if (is.numeric(raw)) as.numeric(raw) else as.numeric(factor(raw))
    limits <- range(numeric_raw, na.rm = TRUE)
    scaled <- if (!all(is.finite(limits)) || diff(limits) == 0) {
      rep(0.5, length(numeric_raw))
    } else {
      (numeric_raw - limits[1]) / diff(limits)
    }
    data.frame(
      feature = feature,
      shap = as.numeric(sv$S[, feature]),
      value = scaled,
      y = i + runif(length(scaled), -0.22, 0.22)
    )
  }))
  p_bee <- ggplot(bee, aes(shap, y, colour = value)) +
    geom_vline(xintercept = 0, colour = "#A9A9A9", linewidth = 0.55) +
    geom_point(size = 1.55, alpha = 0.82, stroke = 0) +
    scale_y_continuous(breaks = seq_along(order_bottom),
                       labels = pretty_label(order_bottom),
                       expand = expansion(add = 0.6)) +
    scale_colour_gradientn(colours = value_colours, limits = c(0, 1),
                           breaks = c(0, 1), labels = c("Low", "High"),
                           name = "Feature value") +
    labs(x = "SHAP value (impact on predicted risk)", y = NULL,
         title = paste0(model_name, " — SHAP summary"),
         subtitle = "Positive values increase predicted risk; negative values decrease it") +
    theme_shap() + theme(legend.position = "right")
  save_publication(p_bee, file.path(model_output, "SHAP_summary"),
                   7.2, max(4.2, 0.28 * nrow(importance) + 2))

  dependence_dir <- file.path(model_output, "dependence")
  dir.create(dependence_dir, recursive = TRUE, showWarnings = FALSE)
  for (feature in colnames(sv$S)) {
    p_dep <- sv_dependence(sv, v = feature, color_var = NULL) +
      geom_hline(yintercept = 0, colour = "#A9A9A9", linewidth = 0.45) +
      labs(title = paste0(model_name, " — ", pretty_label(feature)),
           subtitle = "SHAP dependence in the validation set 2 cohort",
           x = pretty_label(feature), y = "SHAP value") +
      theme_shap()
    safe_name <- substr(gsub("[^[:alnum:]_.-]+", "_", feature), 1, 120)
    save_publication(p_dep, file.path(dependence_dir, paste0(safe_name, "_dependence")),
                     5.4, 4.2)
  }

  local_index_file <- file.path(model_input, "local", "local_explanation_rows.csv")
  if (file.exists(local_index_file)) {
    local_index <- read.csv(local_index_file, check.names = FALSE)
    for (i in seq_len(nrow(local_index))) {
      label <- as.character(local_index$explanation[i])
      row_id <- as.integer(local_index$shap_row[i])
      p_local <- sv_waterfall(sv, row_id = row_id, max_display = ncol(sv$S)) +
        labs(title = paste0(model_name, " — ", gsub("_", " ", label), " patient"),
             subtitle = "Local SHAP explanation of predicted risk") +
        theme_shap()
      save_publication(p_local,
                       file.path(model_output, "local", paste0(label, "_waterfall")),
                       7.2, max(4.2, 0.28 * ncol(sv$S) + 2))
    }
  }
}

comparison <- bind_rows(importance_list) %>%
  mutate(
    model = factor(model, levels = models),
    variable_label = canonical_label(variable),
    variable_label = forcats::fct_reorder(variable_label, importance_share)
  )
p_comparison <- ggplot(comparison, aes(importance_share, variable_label, colour = model)) +
  geom_point(position = position_dodge(width = 0.55), size = 2.3) +
  scale_colour_manual(values = model_colours) +
  labs(x = "Within-model share of mean absolute SHAP", y = NULL, colour = "Model",
       title = "SHAP importance comparison",
       subtitle = "Importance is normalized within each model") +
  theme_shap() + theme(legend.position = "top")
save_publication(p_comparison, file.path(output_root, "SHAP_model_comparison"),
                 7.2, max(4.2, 0.32 * length(unique(comparison$variable_label)) + 2))

write.csv(bind_rows(importance_list),
          file.path(output_root, "SHAP_importance_source_data.csv"), row.names = FALSE)
cat("Publication SHAP figures written to: ", normalizePath(output_root), "\n")


Charcutpoint <- surv_cutpoint(traindata, time = "OS", event = "status_os",
                              variables = "sagittal_max_2d_fd")
fd_score <- as.numeric(Charcutpoint["sagittal_max_2d_fd"][[1]]["estimate"][[1]][[1]])
traindata$min_FD_featureG <- ifelse(traindata$sagittal_max_2d_fd >= fd_score, "High", "Low")
validation_data1$min_FD_featureG <- ifelse(validation_data1$sagittal_max_2d_fd >= fd_score, "High", "Low")
validation_data2_model_data$min_FD_featureG <- ifelse(validation_data2_model_data$sagittal_max_2d_fd >= fd_score, "High", "Low")
survfit(Surv(OS,status_os)~min_FD_featureG, data = traindata) %>%
  ggsurvplot(
    conf.int = TRUE,
    xlab = "Time (months)",
    ylab = "Cumulative Survival Probability",
    main = "Survival Curve for Principal Curvature Min Max",
    legend.title = "Principal Curvature Min Max",
    pval = TRUE,
    palette = c("#00468B", "#ED0000"),
    risk.table = TRUE,
    xlim = c(0,72),
    break.time.by = 12
  )
survfit(Surv(OS,status_os)~min_FD_featureG, data = validation_data1) %>%
  ggsurvplot(
    conf.int = TRUE,
    xlab = "Time (months)",
    ylab = "Cumulative Survival Probability",
    main = "Survival Curve for Principal Curvature Min Max",
    legend.title = "Principal Curvature Min Max",
    pval = TRUE,
    palette = c("#00468B", "#ED0000"),
    risk.table = TRUE,
    xlim = c(0,72),
    break.time.by = 12
  )
survfit(Surv(OS,status_os)~min_FD_featureG, data = validation_data2_model_data) %>%
  ggsurvplot(
    conf.int = TRUE,
    xlab = "Time (months)",
    ylab = "Cumulative Survival Probability",
    main = "Survival Curve for Principal Curvature Min Max",
    legend.title = "Principal Curvature Min Max",
    pval = TRUE,
    palette = c("#00468B", "#ED0000"),
    risk.table = TRUE,
    xlim = c(0,72),
    break.time.by = 12
  )

coxph(Surv(OS, status_os) ~sagittal_max_2d_fd, data = traindata) %>% summary()
coxph(Surv(OS, status_os) ~min_FD_featureG, data = traindata) %>% summary()
traindata$sagittal_max_2d_fd_z <- scale(traindata$sagittal_max_2d_fd)[, 1]

mean(traindata$sagittal_max_2d_fd_z)
sd(traindata$sagittal_max_2d_fd_z)

fit_z <- coxph(Surv(OS, status_os) ~ sagittal_max_2d_fd_z, data = traindata)
summary(fit_z)

# 校准曲线补充指标：Brier 评分、校准斜率、校准截距
# 指标定义：
#   1) Brier score            预测概率与观测结局的均方误差，越小越好；
#   2) Calibration slope      以 glm(y ~ logit(p)) 估计，= 1 表示校准良好，
#                             < 1 提示过拟合 / 风险过度分散，> 1 提示过度集中；
#   3) Calibration intercept  以 glm(y ~ 1 + offset(logit(p))) 估计，
#                             = 0 表示校准良好（calibration-in-the-large）；
#   4) Brier_IPCW             以 survival::brier() 计算的 IPCW 加权 Brier 评分（敏感性指标）。
#
# 评价时点与脚本中已有的校准曲线一致：12 / 24 / 36 个月。
# 每个时点剔除“该时点之前删失、结局未知”的个体，仅用结局确定的个体构造二分类终点：
#   y = I(OS <= t0 且 status_os == 1)
# 预测事件概率 p = 1 - S(t0 | risk score)，其中风险评分在各队列内用 Cox 模型重新校准
# （与已有校准曲线中 rms::cph(Surv(OS, status_os) ~ risk_xgboost) 的做法一致）。

if (!exists("traindata") || !exists("validation_data2_model_data") || !exists("validation_data1")) {
  stop(
    "未找到 traindata / validation_data2_model_data / validation_data1。\n",
    "请先将 SHAP_ONLY 设为 FALSE 并运行本脚本的主分析部分，再执行校准补充指标的计算。"
  )
}

calibration_dir <- file.path(output_dir, "calibration_brier")
dir.create(calibration_dir, recursive = TRUE, showWarnings = FALSE)

cal_times <- c(12, 24, 36)
cal_boot_B <- 1000

cal_model_labels <- c(
  risk_coxph   = "CoxPH",
  risk_gbm     = "GBM",
  risk_lasso   = "LASSO",
  risk_rpart   = "RPart",
  risk_rsf     = "RSF",
  risk_svm     = "SVM",
  risk_xgboost = "XGBoost"
)
cal_risk_vars <- names(cal_model_labels)
cal_risk_vars <- cal_risk_vars[cal_risk_vars %in% names(traindata)]

# 数值化：兼容 haven_labelled / factor
cal_numeric <- function(x) {
  if (inherits(x, "haven_labelled") && requireNamespace("haven", quietly = TRUE)) {
    x <- haven::zap_labels(x)
  }
  if (is.factor(x)) {
    return(suppressWarnings(as.numeric(as.character(x))))
  }
  as.numeric(x)
}

# 固定时点二分类终点：剔除该时点之前删失（结局未知）的个体
cal_binary_horizon <- function(time, status, horizon) {
  time <- cal_numeric(time)
  status <- cal_numeric(status)
  had_event <- time <= horizon & status == 1
  known <- had_event | time > horizon
  list(
    keep = known & is.finite(time) & is.finite(status),
    y = ifelse(had_event, 1, 0)
  )
}

# 从 summary.survfit 的结果中取回全部个体的生存概率
cal_extract_surv <- function(surv_values, n_subject) {
  if (is.matrix(surv_values)) {
    if (ncol(surv_values) == n_subject) {
      return(as.numeric(surv_values[1, ]))
    }
    if (nrow(surv_values) == n_subject) {
      return(as.numeric(surv_values[, 1]))
    }
  }
  as.numeric(surv_values)
}

# 关注时点的预测事件概率：p = 1 - S(t0 | risk score)
cal_predict_risk <- function(cox_model, newdata, horizon) {
  sf <- survival::survfit(cox_model, newdata = newdata)
  sm <- summary(sf, times = horizon, extend = TRUE)
  surv_hat <- cal_extract_surv(sm$surv, nrow(newdata))
  if (length(surv_hat) != nrow(newdata)) {
    stop("预测生存概率个数(", length(surv_hat), ")与个体数(", nrow(newdata), ")不一致。")
  }
  pmin(pmax(1 - surv_hat, 1e-6), 1 - 1e-6)
}

cal_na_row <- function(cohort_label, model_label, t0, n, events, brier = NA_real_) {
  data.frame(
    cohort = cohort_label,
    model = model_label,
    time = t0,
    n = n,
    events = events,
    Brier = brier,
    Brier_IPCW = NA_real_,
    calibration_intercept = NA_real_,
    intercept_p = NA_real_,
    calibration_slope = NA_real_,
    slope_lower_95CI = NA_real_,
    slope_upper_95CI = NA_real_,
    slope_p_for_1 = NA_real_,
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}

calibration_brier_table <- function(data, cohort_label,
                                   risk_vars = cal_risk_vars,
                                   horizons = cal_times,
                                   time_var = "OS", status_var = "status_os") {
  time_all <- cal_numeric(data[[time_var]])
  status_all <- cal_numeric(data[[status_var]])

  purrr::map_dfr(risk_vars, function(risk_var) {
    model_label <- unname(cal_model_labels[risk_var])
    risk_all <- as.numeric(data[[risk_var]])
    usable <- !is.na(time_all) & !is.na(status_all) & !is.na(risk_all) &
      is.finite(time_all) & is.finite(status_all) & is.finite(risk_all) & time_all > 0

    cox_dat <- data.frame(
      time = time_all[usable],
      status = status_all[usable],
      score = risk_all[usable]
    )

    if (nrow(cox_dat) < 20 || length(unique(cox_dat$score)) < 2 ||
        sum(cox_dat$status == 1) < 5) {
      return(purrr::map_dfr(horizons, function(t0) {
        cal_na_row(cohort_label, model_label, t0, nrow(cox_dat), sum(cox_dat$status == 1))
      }))
    }

    cox_fit <- survival::coxph(
      survival::Surv(time, status) ~ score,
      data = cox_dat, x = TRUE, y = TRUE
    )

    purrr::map_dfr(horizons, function(t0) {
      bin <- cal_binary_horizon(time_all, status_all, t0)
      keep <- bin$keep & usable
      keep[is.na(keep)] <- FALSE
      y <- bin$y[keep]

      risk_pred <- tryCatch(
        cal_predict_risk(cox_fit, data.frame(score = risk_all[keep]), t0),
        error = function(e) NULL
      )
      if (is.null(risk_pred) || length(unique(y)) < 2 ||
          length(unique(round(risk_pred, 8))) < 2) {
        return(cal_na_row(
          cohort_label, model_label, t0, length(y), sum(y == 1),
          brier = if (is.null(risk_pred)) NA_real_ else mean((y - risk_pred)^2)
        ))
      }

      brier_ipcw <- tryCatch({
        b <- survival::brier(cox_fit, times = t0, newdata = cox_dat)
        v <- as.numeric(b$brier)
        v[length(v)]
      }, error = function(e) NA_real_)

      lp <- stats::qlogis(risk_pred)
      intercept_fit <- stats::glm(y ~ 1 + offset(lp), family = stats::binomial())
      slope_fit <- stats::glm(y ~ lp, family = stats::binomial())
      slope_ci <- suppressMessages(stats::confint.default(slope_fit, parm = "lp"))
      slope_est <- stats::coef(slope_fit)[["lp"]]
      slope_se <- summary(slope_fit)$coefficients["lp", "Std. Error"]
      slope_p_for_1 <- if (is.finite(slope_se) && slope_se > 0) {
        2 * stats::pnorm(-abs((slope_est - 1) / slope_se))
      } else {
        NA_real_
      }

      data.frame(
        cohort = cohort_label,
        model = model_label,
        time = t0,
        n = length(y),
        events = sum(y == 1),
        Brier = mean((y - risk_pred)^2),
        Brier_IPCW = brier_ipcw,
        calibration_intercept = unname(stats::coef(intercept_fit)[1]),
        intercept_p = unname(summary(intercept_fit)$coefficients[1, "Pr(>|z|)"]),
        calibration_slope = unname(slope_est),
        slope_lower_95CI = unname(slope_ci[1]),
        slope_upper_95CI = unname(slope_ci[2]),
        slope_p_for_1 = unname(slope_p_for_1),
        row.names = NULL,
        stringsAsFactors = FALSE
      )
    })
  })
}

calibration_results <- dplyr::bind_rows(
  calibration_brier_table(traindata, "Training set"),
  calibration_brier_table(validation_data2_model_data, "Validation set 2"),
  calibration_brier_table(validation_data1, "Validation set 1")
)

write.csv(
  calibration_results,
  file.path(calibration_dir, "calibration_Brier_slope_intercept_all_models.csv"),
  row.names = FALSE
)

# 整理为可直接入稿的汇总表
calibration_summary <- calibration_results %>%
  dplyr::mutate(
    `Time (months)` = time,
    `Brier score` = round(Brier, 3),
    `Brier score (IPCW)` = round(Brier_IPCW, 3),
    `Calibration intercept` = round(calibration_intercept, 3),
    `Intercept P` = signif(intercept_p, 3),
    `Calibration slope (95% CI)` = ifelse(
      is.na(calibration_slope), NA_character_,
      sprintf("%.3f (%.3f-%.3f)", calibration_slope, slope_lower_95CI, slope_upper_95CI)
    ),
    `Slope P for 1` = signif(slope_p_for_1, 3)
  ) %>%
  dplyr::select(
    cohort, model, `Time (months)`, n, events,
    `Brier score`, `Brier score (IPCW)`,
    `Calibration intercept`, `Intercept P`,
    `Calibration slope (95% CI)`, `Slope P for 1`
  ) %>%
  dplyr::arrange(cohort, model, `Time (months)`)

write.csv(
  calibration_summary,
  file.path(calibration_dir, "calibration_summary_table_all_models.csv"),
  row.names = FALSE
)

# 校准曲线所对应的模型（XGBoost）单独出一张表
calibration_summary_xgboost <- calibration_summary %>%
  dplyr::filter(model == "XGBoost") %>%
  dplyr::select(-model)

write.csv(
  calibration_summary_xgboost,
  file.path(calibration_dir, "calibration_summary_table_XGBoost.csv"),
  row.names = FALSE
)

print(calibration_summary_xgboost)

# 将 Brier 评分、校准斜率、校准截距直接标注到校准曲线上（XGBoost，12/24/36 个月）
annotate_calibration_curves <- function(data, cohort_label, cohort_tag,
                                        risk_var = "risk_xgboost",
                                        horizons = cal_times,
                                        m = NULL, B = cal_boot_B,
                                        out_dir = calibration_dir) {
  cal_dat <- data.frame(
    cal_time = cal_numeric(data[["OS"]]),
    cal_status = cal_numeric(data[["status_os"]]),
    risk_score = as.numeric(data[[risk_var]])
  )
  cal_dat <- cal_dat[stats::complete.cases(cal_dat), , drop = FALSE]

  if (nrow(cal_dat) < 20 || length(unique(cal_dat$risk_score)) < 2 ||
      sum(cal_dat$cal_status == 1) < 5) {
    warning(cohort_label, ": 可用于校准分析的数据不足，跳过标注图。")
    return(invisible(NULL))
  }
  if (is.null(m)) m <- max(20L, min(60L, floor(nrow(cal_dat) / 4)))

  old_dd <- getOption("datadist")
  options(datadist = rms::datadist(cal_dat))
  on.exit(options(datadist = old_dd), add = TRUE)

  cols <- c("#3B83C0", "#EE362B", "#F8CB75")
  curves <- vector("list", length(horizons))
  names(curves) <- as.character(horizons)

  for (i in seq_along(horizons)) {
    t0 <- horizons[[i]]
    fit <- rms::cph(
      survival::Surv(cal_time, cal_status) ~ risk_score,
      data = cal_dat, x = TRUE, y = TRUE, surv = TRUE, time.inc = t0
    )
    curves[[i]] <- tryCatch(
      rms::calibrate(fit, cmethod = "KM", method = "boot", u = t0, m = m, B = B),
      error = function(e) {
        auto_m <- max(20L, min(60L, floor(nrow(cal_dat) / 4)))
        warning("m = ", m, " 时 calibrate 失败（", conditionMessage(e),
                "），改用 m = ", auto_m, "。")
        rms::calibrate(fit, cmethod = "KM", method = "boot", u = t0, m = auto_m, B = B)
      }
    )
  }

  ann_txt <- paste(
    c(
      cohort_label,
      vapply(horizons, function(t0) {
        row <- calibration_results[
          calibration_results$cohort == cohort_label &
            calibration_results$model == "XGBoost" &
            calibration_results$time == t0, , drop = FALSE
        ]
        if (nrow(row) == 0 || is.na(row$calibration_slope[1])) {
          return(sprintf("%dm: Brier = NA, intercept = NA, slope = NA", t0))
        }
        sprintf(
          "%dm: Brier = %.3f, intercept = %.2f, slope = %.2f (%.2f-%.2f)",
          t0, row$Brier[1], row$calibration_intercept[1],
          row$calibration_slope[1], row$slope_lower_95CI[1], row$slope_upper_95CI[1]
        )
      }, character(1))
    ),
    collapse = "\n"
  )

  draw <- function() {
    plot(curves[[1]], lwd = 2, lty = 1, errbar.col = cols[1],
         xlim = c(0, 1), ylim = c(0, 1),
         xlab = "Predicted OS probability", ylab = "Observed OS probability",
         col = cols[1])
    if (length(curves) > 1) {
      for (i in 2:length(curves)) {
        plot(curves[[i]], add = TRUE, lwd = 2, lty = 1,
             errbar.col = cols[i], col = cols[i])
      }
    }
    legend("bottomright", paste0(horizons / 12, "-year OS"),
           lwd = 2, col = cols[seq_along(horizons)], bty = "n", cex = 0.85)
    graphics::text(x = 0.02, y = 0.98, labels = ann_txt, adj = c(0, 1), cex = 0.72)
  }

  grDevices::pdf(file.path(out_dir, paste0("Calibration_annotated_", cohort_tag, ".pdf")), 7, 7)
  draw()
  grDevices::dev.off()

  grDevices::tiff(
    file.path(out_dir, paste0("Calibration_annotated_", cohort_tag, ".tif")),
    width = 7, height = 7, units = "in", res = 300, compression = "lzw"
  )
  draw()
  grDevices::dev.off()

  invisible(curves)
}

annotated_calibration_plots <- list(
  "Training set" = tryCatch(
    annotate_calibration_curves(traindata, "Training set", "TrainingSet", m = 60),
    error = function(e) {
      warning("训练集校准标注图失败：", conditionMessage(e))
      NULL
    }
  ),
  "Validation set 2" = tryCatch(
    annotate_calibration_curves(validation_data2_model_data, "Validation set 2", "ValidationSet2", m = 27),
    error = function(e) {
      warning("验证集2校准标注图失败：", conditionMessage(e))
      NULL
    }
  ),
  "Validation set 1" = tryCatch(
    annotate_calibration_curves(
      validation_data1, "Validation set 1", "ValidationSet1", m = 27
    ),
    error = function(e) {
      warning("验证集1校准标注图失败：", conditionMessage(e))
      NULL
    }
  )
)

save(
  calibration_results,
  calibration_summary,
  calibration_summary_xgboost,
  annotated_calibration_plots,
  file = file.path(calibration_dir, "calibration_Brier_slope_intercept.RData")
)

cat(
  "\n校准补充指标（Brier 评分、校准斜率、校准截距）已保存至：\n",
  normalizePath(calibration_dir, mustWork = FALSE), "\n", sep = ""
)
# Biomarker interpretability
Clinical_Data <- read_sav(development_clinical_file)
FD_Data <- read.csv(development_fd_file)
FD_Data$PatientsNum <- as.numeric(str_sub(FD_Data$patient_id,start = 0, end = 7))
merged_data <- merge(Clinical_Data, FD_Data[,-which(colnames(FD_Data) == "patient_id")], by = "PatientsNum")
merged_data$BM <- ifelse(merged_data$BM == 1, 1, 0)
merged_data$status_os <- ifelse(merged_data$status_os == 1, 1, 0)
merged_data[,c(2,4,12,15,39,41:42,45)] <- lapply(merged_data[,c(2,4,12,15,39,41:42,45)],function(x){
  factor(x,levels = c(data.frame(Hmisc::valueLabel(x))[,1]),labels = rownames(data.frame(Hmisc::valueLabel(x))))
})
non_pcidata <- merged_data
non_pcidata$ECOG <- ifelse(non_pcidata$ECOG>1,">1","0-1")
non_pcidata$PCI <- ifelse(non_pcidata$PCI==0,"No", "Yes")
non_pcidata$RT <- ifelse(non_pcidata$TimeToRT==0,"No", "Yes")
non_pcidata$Chemotherapy <- ifelse(non_pcidata$Chemotherapy>=4,"Yes", "No")

non_pcidata$PatientsNum_chr <- as.character(non_pcidata$PatientsNum)
non_pcidata$year <- as.integer(substr(non_pcidata$PatientsNum_chr, start = 1, stop = 2))+2000
biomarkers <- c(
  "LDH",  "NSE", "Cyfra211", "ProGRP"
)
MoLan_Data1 <- non_pcidata[,c(names(non_pcidata[c(1,2:4,8,15,40,50,39,67:82,64,63,85)]),biomarkers)]
MoLan_Data1 <- MoLan_Data1[!MoLan_Data1$TreatmentMode %in% c("CT_IM" ,"CT_IM_RT"),]
MoLan_Data1$TreatmentMode <- ifelse(MoLan_Data1$TreatmentMode %in% c("CRT","C_CRT_C"),"CRT",
                                    ifelse(MoLan_Data1$TreatmentMode=="C","C","CTRT"))
MoLan_Data1$group <- "shanxi"


FD_Data_tianjin <- read.csv(external_fd_file)
FD_Data_tianjin$PatientsNum <- as.numeric(str_extract(FD_Data_tianjin$patient_id, "^\\d+"))
clinical_data_tianjin <- read_xlsx(external_clinical_file)
clinical_data_tianjin$PatientsNum <- clinical_data_tianjin$序号
clinical_data_tianjin$ECOG <- ifelse(clinical_data_tianjin$`KPS≥80`<=1,"0-1",">1")
clinical_data_tianjin$Chemotherapy <- ifelse(clinical_data_tianjin$化疗总周期数<4,"No","Yes")
clinical_data_tianjin$RT <- "Yes"
clinical_data_tianjin$PCI <- ifelse(clinical_data_tianjin$PCI==1,"Yes","No")
clinical_data_tianjin <- merge.data.frame(clinical_data_tianjin,FD_Data_tianjin,by = "PatientsNum")

clinical_data_tianjin$status_os <- clinical_data_tianjin$状态
clinical_data_tianjin <- clinical_data_tianjin[!(clinical_data_tianjin$OS<6&clinical_data_tianjin$status_os==0),]
clinical_data_tianjin$Age <- clinical_data_tianjin$年龄
clinical_data_tianjin$Gender <- ifelse(clinical_data_tianjin$性别==1,"Male","Female")
clinical_data_tianjin$Clinical_Stage <- ifelse(clinical_data_tianjin$`临床分期（1:I期;2:II期;3:IIIA期;4:IIIB期;5:IIIC期）`==3,"IIIa",
                                               ifelse(clinical_data_tianjin$`临床分期（1:I期;2:II期;3:IIIA期;4:IIIB期;5:IIIC期）`==4,"IIIb",
                                                      ifelse(clinical_data_tianjin$`临床分期（1:I期;2:II期;3:IIIA期;4:IIIB期;5:IIIC期）`==5,"III‎c","I-IIb")))
clinical_data_tianjin$Smoke <- ifelse(clinical_data_tianjin$吸烟==1,"Yes","No")
clinical_data_tianjin$TreatmentMode <- ifelse(clinical_data_tianjin$同步化疗有无==1, "CRT","CTRT")
clinical_data_tianjin$ITHScore <- clinical_data_tianjin$ith_score
clinical_data_tianjin$date_correct_base <- as.Date(clinical_data_tianjin$确诊时间, format = "%m/%d/%Y")
clinical_data_tianjin$year <- as.integer(format(clinical_data_tianjin$date_correct_base, "%Y"))
clinical_data_tianjin$group <- "tianjin"

clinical_data_tianjin <- clinical_data_tianjin %>%
  mutate(
    NEU      = 治疗前中性粒细胞,
    LYM      = 治疗前淋巴细胞,
    HGB      = 治疗前血红蛋白,
    Alb      = 治疗前ALB,
    LDH      = 治疗前LDH,
    CEA      = 治疗前CEA,
    CA199    = 治疗前CA199,
    NSE      = 治疗前NSE,
    Cyfra211 = 治疗前Cyfra21,
    ProGRP   = 治疗前ProGRP
  )

MoLan_Data2 <- clinical_data_tianjin[names(MoLan_Data1)]
MoLan_Data <- rbind.data.frame(MoLan_Data1,MoLan_Data2)
MoLan_Data$group <- ifelse(MoLan_Data$year>2018,">2018","2013-2018")
MoLan_Data <- na.omit(MoLan_Data)
MoLan_Data$PatientsNum <- paste0(MoLan_Data$OS,MoLan_Data$sagittal_max_2d_fd)
traindata$PatientsNum <- paste0(traindata$OS,traindata$sagittal_max_2d_fd)
validation_data1$PatientsNum <- paste0(validation_data1$OS,validation_data1$sagittal_max_2d_fd)
validation_data2_model_data$PatientsNum <- paste0(validation_data2_model_data$OS,validation_data2_model_data$sagittal_max_2d_fd)
traindata$cohort <- "Training Cohort"
traindata <- merge.data.frame(traindata,MoLan_Data,by = "PatientsNum")
validation_data1$cohort <- "Validation Cohort 1"
validation_data1 <- merge.data.frame(validation_data1,MoLan_Data,by = "PatientsNum")
validation_data2_model_data$cohort <- "Validation Cohort 2"
validation_data2_model_data <- merge.data.frame(validation_data2_model_data,MoLan_Data,by = "PatientsNum")
cohort_data <- bind_rows(traindata[,c(biomarkers,"cohort","xgboost_group","risk_xgboost","PatientsNum")],
                         validation_data1[,c(biomarkers,"cohort","xgboost_group","risk_xgboost","PatientsNum")],
                         validation_data2_model_data[,c(biomarkers,"cohort","xgboost_group","risk_xgboost","PatientsNum")])
cohort_data$xgboost_group <- as.factor(cohort_data$xgboost_group)
# Combined is descriptive only; inferential validation should use each cohort separately.
combined_data <- cohort_data %>% mutate(cohort = "Combined (descriptive)")

combined_data1 <- cohort_data[cohort_data$cohort != "Training Cohort",] %>%
  mutate(cohort = "Validation sets 1 + 2")

analysis_data <- bind_rows(cohort_data, combined_data)

analysis_data <- bind_rows(analysis_data, combined_data1)

write.csv(cohort_data, file.path(output_dir, "cohort_data_joined_by_PatientsNum.csv"),
          row.names = FALSE)
write.csv(
  cohort_data %>% count(cohort, name = "n_patients"),
  file.path(output_dir, "cohort_sample_sizes.csv"), row.names = FALSE
)

safe_numeric <- function(x) suppressWarnings(as.numeric(as.character(x)))

correlation_one <- function(data, biomarker, method) {
  x <- safe_numeric(data[[biomarker]])
  y <- safe_numeric(data$risk_xgboost)
  ok <- is.finite(x) & is.finite(y)
  n <- sum(ok)
  if (n < 4L || length(unique(x[ok])) < 2L || length(unique(y[ok])) < 2L) {
    return(data.frame(n = n, estimate = NA_real_, conf_low = NA_real_,
                      conf_high = NA_real_, p_value = NA_real_))
  }
  fit <- suppressWarnings(cor.test(x[ok], y[ok], method = method, exact = FALSE))
  ci <- if (!is.null(fit$conf.int)) as.numeric(fit$conf.int) else c(NA_real_, NA_real_)
  data.frame(n = n, estimate = unname(fit$estimate), conf_low = ci[1],
             conf_high = ci[2], p_value = fit$p.value)
}

correlation_results <- bind_rows(lapply(unique(analysis_data$cohort), function(cohort_name) {
  dat <- analysis_data %>% filter(cohort == cohort_name)
  bind_rows(lapply(c("pearson", "spearman"), function(method) {
    bind_rows(lapply(biomarkers, function(marker) {
      cbind(data.frame(cohort = cohort_name, biomarker = marker, method = method),
            correlation_one(dat, marker, method))
    }))
  }))
})) %>%
  group_by(cohort, method) %>%
  mutate(p_fdr_bh = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  mutate(significant_fdr = !is.na(p_fdr_bh) & p_fdr_bh < 0.05)

welch_one <- function(data, biomarker) {
  value <- safe_numeric(data[[biomarker]])
  group <- droplevels(data$xgboost_group)
  ok <- is.finite(value) & !is.na(group)
  value <- value[ok]; group <- droplevels(group[ok])
  n_low <- sum(group == "Low"); n_high <- sum(group == "High")
  if (length(levels(group)) != 2L || n_low < 2L || n_high < 2L) {
    return(data.frame(n_low=n_low, n_high=n_high, mean_low=NA, mean_high=NA,
                      sd_low=NA, sd_high=NA, mean_difference=NA, conf_low=NA,
                      conf_high=NA, t_statistic=NA, df=NA, p_value=NA,
                      hedges_g=NA))
  }
  low <- value[group == "Low"]; high <- value[group == "High"]
  fit <- t.test(high, low, var.equal = FALSE)
  pooled_sd <- sqrt(((n_high - 1) * var(high) + (n_low - 1) * var(low)) /
                      (n_high + n_low - 2))
  cohen_d <- if (is.finite(pooled_sd) && pooled_sd > 0) (mean(high) - mean(low)) / pooled_sd else NA_real_
  correction <- 1 - 3 / (4 * (n_high + n_low) - 9)
  data.frame(
    n_low=n_low, n_high=n_high, mean_low=mean(low), mean_high=mean(high),
    sd_low=sd(low), sd_high=sd(high), mean_difference=mean(high)-mean(low),
    conf_low=fit$conf.int[1], conf_high=fit$conf.int[2],
    t_statistic=unname(fit$statistic), df=unname(fit$parameter),
    p_value=fit$p.value, hedges_g=cohen_d * correction
  )
}
ttest_results <- bind_rows(lapply(unique(analysis_data$cohort), function(cohort_name) {
  dat <- analysis_data %>% filter(cohort == cohort_name)
  bind_rows(lapply(biomarkers, function(marker) {
    cbind(data.frame(cohort = cohort_name, biomarker = marker), welch_one(dat, marker))
  }))
})) %>%
  group_by(cohort) %>%
  mutate(p_fdr_bh = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  mutate(significant_fdr = !is.na(p_value) & p_value < 0.05)
ttest_star <- ttest_results %>%
  filter(cohort %in% c("Training Cohort", "Validation Cohort 1", "Validation Cohort 2", "Validation sets 1 + 2")) %>%
  mutate(
    biomarker = factor(biomarker, levels = biomarkers),
    star = ifelse(significant_fdr, "*", "")
  )

# Biomarker boxplots with FDR significance markers.
theme_bio <- theme_classic(base_size = 8, base_family = "Arial") +
  theme(axis.text = element_text(colour = "black"),
        strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"))
long_box <- analysis_data %>%
  dplyr::select(cohort, PatientsNum, xgboost_group, dplyr::all_of(biomarkers)) %>%
  pivot_longer(dplyr::all_of(biomarkers), names_to = "biomarker", values_to = "value") %>%
  mutate(
    value = safe_numeric(value),
    biomarker = factor(biomarker, levels = biomarkers)
  )

# 合并星号数据
long_box <- left_join(long_box, ttest_star, by = c("cohort", "biomarker"))

# 计算标注y位置：取每组最大值+一点偏移，避免盖住离群点
star_anno <- long_box %>%
  group_by(cohort, biomarker) %>%
  summarise(y_max = max(value, na.rm = TRUE) * 1.03, star = first(star), .groups = "drop")

p_box <- ggplot(long_box[long_box$cohort %in% c("Training Cohort", "Validation sets 1 + 2"),], aes(x = xgboost_group, y = value, fill = xgboost_group)) +
  geom_boxplot(width = .62, outlier.shape = NA, linewidth = .3) +
  geom_jitter(width = .14, size = .35, alpha = .35) +
  # 重点：facet_grid，cohort分行，biomarker分列，替代原来facet_wrap
  facet_grid(
    rows = vars(cohort),
    cols = vars(biomarker),
    scales = "free_y"
  ) +
  # 添加显著性*
  geom_text(data = star_anno[star_anno$cohort %in% c("Training Cohort", "Validation sets 1 + 2"),], aes(x = 1.5, y = y_max, label = star),
            inherit.aes = FALSE, size = 2.8, fontface = "bold", color = "#C94235") +
  scale_fill_manual(values = c(Low = "#3568A8", High = "#C94235")) +
  labs(
    title = "Laboratory biomarkers by ITHScore group",
    subtitle = "Each cohort–biomarker panel has an independent y-axis; * BH-FDR < 0.05",
    x = NULL, y = NULL
  ) + theme_bio +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle=35, hjust=1),
    strip.text = element_text(size = 6.2)
  )
p_box
ggplot(long_box[long_box$cohort %in% c("Training Cohort", "Validation sets 1 + 2") & long_box$biomarker %in% c("ProGRP"),], aes(x = xgboost_group, y = value, fill = xgboost_group)) +
  geom_boxplot(width = .62, outlier.shape = NA, linewidth = .3) +
  geom_jitter(width = .14, size = .35, alpha = .35) +
  # 重点：facet_grid，cohort分行，biomarker分列，替代原来facet_wrap
  facet_grid(
    rows = vars(cohort),
    cols = vars(biomarker),
    scales = "free_y"
  ) +
  # 添加显著性*
  geom_text(data = star_anno[star_anno$cohort %in% c("Training Cohort", "Validation sets 1 + 2") & star_anno$biomarker %in% c("ProGRP"),], aes(x = 1.5, y = y_max, label = star),
            inherit.aes = FALSE, size = 2.8, fontface = "bold", color = "#C94235") +
  scale_fill_manual(values = c(Low = "#3568A8", High = "#C94235")) +
  labs(
    title = "Laboratory biomarkers by ITHScore group",
    subtitle = "Each cohort–biomarker panel has an independent y-axis; * BH-FDR < 0.05",
    x = NULL, y = NULL
  ) + theme_bio +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle=35, hjust=1),
    strip.text = element_text(size = 6.2)
  )
# Standardized biomarker boxplots.
long_box_z <- long_box[long_box$cohort %in% c("Training Cohort", "Validation sets 1 + 2"),] %>%
  group_by(cohort, biomarker) %>%
  mutate(
    value_z = if (sum(is.finite(value)) >= 2L && stats::sd(value, na.rm = TRUE) > 0) {
      as.numeric(scale(value))
    } else {
      NA_real_
    }
  ) %>%
  ungroup()
p_box_z <- ggplot(long_box_z[long_box_z$cohort %in% c("Training Cohort", "Validation sets 1 + 2"),], aes(biomarker, value_z, fill = xgboost_group)) +
  geom_hline(yintercept = 0, colour = "#A0A0A0", linewidth = .35) +
  geom_boxplot(
    position = position_dodge(width = .72), width = .62,
    outlier.shape = NA, linewidth = .3
  ) +
  # 星号标注：放在箱线组上方
  geom_text(
    data = ttest_star[ttest_star$cohort %in% c("Training Cohort", "Validation sets 1 + 2"),],
    aes(x = biomarker, y = 2.6, label = star),
    inherit.aes = FALSE, size = 3, fontface = "bold", color="#C94235"
  ) +
  facet_wrap(~cohort, ncol = 1) +
  scale_fill_manual(values = c(Low = "#3568A8", High = "#C94235")) +
  scale_y_continuous(limits = c(-3, 3), breaks = seq(-3,3,1)) + # 新增y轴范围控制
  labs(
    title = "Standardized laboratory biomarkers by ITHScore group",
    subtitle = "Within-cohort biomarker z-scores; intended for cross-marker visual comparison; * BH-FDR < 0.05",
    x = NULL, y = "Standardized value (z-score)", fill = "ITHScore group"
  ) +
  theme_bio +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 55, hjust = 1, size = 6)
  )

p_box_z
