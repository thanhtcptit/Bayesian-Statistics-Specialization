*WIP*

# Syllabus

## Course 1: From Concept to Data Analysis

This course introduces the Bayesian approach to statistics, starting with the concept of probability and moving to the analysis of data. We will learn about the philosophy of the Bayesian approach as well as how to implement it for common types of data. We will compare the Bayesian approach to the more commonly-taught Frequentist approach, and see some of the benefits of the Bayesian approach. In particular, the Bayesian approach allows for better accounting of uncertainty, results that have more intuitive and interpretable meaning, and more explicit statements of assumptions. This course combines lecture videos, computer demonstrations, readings, exercises, and discussion boards to create an active learning experience. For computing, you have the choice of using Microsoft Excel or the open-source, freely available statistical package R, with equivalent content for both options. The lectures provide some of the basic mathematical development as well as explanations of philosophy and interpretation. Completion of this course will give you an understanding of the concepts of the Bayesian approach, understanding the key differences between Bayesian and Frequentist approaches, and the ability to do basic data analyses.

**Week 1: Probability and Bayes' Theorem**

In this module, we review the basics of probability and Bayes’ theorem. In Lesson 1, we introduce the different paradigms or definitions of probability and discuss why probability provides a coherent framework for dealing with uncertainty. In Lesson 2, we review the rules of conditional probability and introduce Bayes’ theorem. Lesson 3 reviews common probability distributions for discrete and continuous random variables.

**Week 2: Statistical Inference**

This module introduces concepts of statistical inference from both frequentist and Bayesian perspectives. Lesson 4 takes the frequentist view, demonstrating maximum likelihood estimation and confidence intervals for binomial data. Lesson 5 introduces the fundamentals of Bayesian inference. Beginning with a binomial likelihood and prior probabilities for simple hypotheses, you will learn how to use Bayes’ theorem to update the prior with data to obtain posterior probabilities. This framework is extended with the continuous version of Bayes theorem to estimate continuous model parameters, and calculate posterior probabilities and credible intervals.

**Week 3: Priors and Models for Discrete Data**

In this module, you will learn methods for selecting prior distributions and building models for discrete data. Lesson 6 introduces prior selection and predictive distributions as a means of evaluating priors. Lesson 7 demonstrates Bayesian analysis of Bernoulli data and introduces the computationally convenient concept of conjugate priors. Lesson 8 builds a conjugate model for Poisson data and discusses strategies for selection of prior hyperparameters.

**Week 4: Models for Continuous Data**

This module covers conjugate and objective Bayesian analysis for continuous data. Lesson 9 presents the conjugate model for exponentially distributed data. Lesson 10 discusses models for normally distributed data, which play a central role in statistics. In Lesson 11, we return to prior selection and discuss ‘objective’ or ‘non-informative’ priors. Lesson 12 presents Bayesian linear regression with non-informative priors, which yield results comparable to those of classical regression.


## Course 2: Techniques and Models

This is the second of a two-course sequence introducing the fundamentals of Bayesian statistics. It builds on the course Bayesian Statistics: From Concept to Data Analysis, which introduces Bayesian methods through use of simple conjugate models. Real-world data often require more sophisticated models to reach realistic conclusions. This course aims to expand our “Bayesian toolbox” with more general models, and computational techniques to fit them. In particular, we will introduce Markov chain Monte Carlo (MCMC) methods, which allow sampling from posterior distributions that have no analytical solution. We will use the open-source, freely available software R (some experience is assumed, e.g., completing the previous course in R) and JAGS (no experience required). We will learn how to construct, fit, assess, and compare Bayesian statistical models to answer scientific questions involving continuous, binary, and count data. This course combines lecture videos, computer demonstrations, readings, exercises, and discussion boards to create an active learning experience. The lectures provide some of the basic mathematical development, explanations of the statistical modeling process, and a few basic modeling techniques commonly used by statisticians. Computer demonstrations provide concrete, practical walkthroughs. Completion of this course will give you access to a wide range of Bayesian analytical tools, customizable to your data.

**Week 1: Statistical modeling and Monte Carlo estimation**

Statistical modeling, Bayesian modeling, Monte Carlo estimation

Learning Objectives
- Understand how to interpret and specify the components of Bayesian statistical models: likelihood, prior, posterior.
- Given a scientific problem and data, identify an appropriate likelihood/model and relevant hypotheses.
- Use R to generate Monte Carlo samples from posterior distributions to obtain: point estimates, interval estimates, posterior probabilities of hypotheses.

**Week 2: Markov chain Monte Carlo (MCMC)**

Metropolis-Hastings, Gibbs sampling, assessing convergence

Learning Objectives
- Understand the basics of MCMC including Metropolis-Hastings and Gibbs sampling.
- Write a statistical model in JAGS and produce posterior samples by calling JAGS from R.
- Assess MCMC output to determine if it is suitable for inference.

**Week 3: Common statistical models**

Linear regression, ANOVA, logistic regression, multiple factor ANOVA

Learning Objectives
- Assess the performance of a statistical model and compare competing models.
- Specify, fit, assess, and use common statistical models for continuous and binary data.
- Use statistical modeling results to draw scientific conclusions.

**Week 4: Count data and hierarchical modeling**

Poisson regression, hierarchical modeling

Learning Objectives
- Specify, fit, assess, and use common statistical models for count data.
- Extend basic statistical models to account for correlated observations using hierarchical models.
- Use statistical modeling results to draw scientific conclusions.

**Week 5: Capstone project**

Peer-reviewed data analysis project

Learning Objectives
- Demonstrate proficiency in all previous learning objectives.
- Efficiently and effectively communicate results of a data analysis.


## Course 3: Mixture Models

Bayesian Statistics: Mixture Models introduces you to an important class of statistical models. The course is organized in five modules, each of which contains lecture videos, short quizzes, background reading, discussion prompts, and one or more peer-reviewed assignments. Statistics is best learned by doing it, not just watching a video, so the course is structured to help you learn through application. 

Some exercises require the use of R, a freely-available statistical software package. A brief tutorial is provided, but we encourage you to take advantage of the many other resources online for learning R if you are interested.

This is an intermediate-level course, and it was designed to be the third in UC Santa Cruz's series on Bayesian statistics, after Herbie Lee's "Bayesian Statistics: From Concept to Data Analysis" and Matthew Heiner's "Bayesian Statistics: Techniques and Models." To succeed in the course, you should have some knowledge of and comfort with calculus-based probability, principles of maximum-likelihood estimation, and Bayesian estimation.

**Week 1: Basic concepts on Mixture Models**

This module defines mixture models, discusses its properties, and develops the likelihood function for a random sample from a mixture model that will be the basis for statistical learning.

Learning Objectives
- List the goals of the course
- Master the basics of the R computing environment that will be used for this course
- Define mixture models
- Compute the expectation and variance of a mixture distribution
- Write down the likelihood function associated with a random sample from a mixture distribution
- Simulate random samples from a mixture distribution

**Week 2: Maximum likelihood estimation for Mixture Models**

Learning Objectives
- List the basic principles behind the EM algorithm for fitting a mixture model
- Derive the EM algorithm for fitting a mixture model
- Implement in the R environment the EM algorithm for fitting a mixture model

**Week 3: Bayesian estimation for Mixture Models**

Learning Objectives
- List the basic principles behind Markov chain Monte Carlo algorithms for fitting mixture models
- Derive Markov chain Monte Carlo algorithms for fitting mixture models
- Implement in the R environment Markov chain Monte Carlo algorithms for fitting mixture models

**Week 4: Applications of Mixture Models**

Learning Objectives
- Use Mixture Models to provide density estimates
- Use Mixture Models to solve clustering (unsupervised classification) problems
- Use Mixture Models to solve (supervised) classification problems

**Week 5: Practical considerations**

Learning Objectives
- Address computational issues associated with algorithms used to fit Mixture Models
- Explain the implications of multimodality issues in the algorithms used to fit Mixture Models
- Estimate the number of components in a Mixture Model


## Course 4: Time Series Analysis

This course for practicing and aspiring data scientists and statisticians. It is the fourth of a four-course sequence introducing the fundamentals of Bayesian statistics. It builds on the course Bayesian Statistics: From Concept to Data Analysis, Techniques and Models, and Mixture models. 

Time series analysis is concerned with modeling the dependency among elements of a sequence of temporally related variables. To succeed in this course, you should be familiar with calculus-based probability, the principles of maximum likelihood estimation, and Bayesian inference. You will learn how to build models that can describe temporal dependencies and how to perform Bayesian inference and forecasting for the models. You will apply what you've learned with the open-source, freely available software R with sample databases. Your instructor Raquel Prado will take you from basic concepts for modeling temporally dependent data to implementation of specific classes of models

**Week 1: Introduction to time series and the AR(1) process**

This module defines stationary time series processes, the autocorrelation function and the autoregressive process of order one or AR(1). Parameter estimation via maximum likelihood and Bayesian inference in the AR(1) are also discussed.

Learning Objectives
- List the goals of the course and identify the basics of the R environment
- Explain stationary time series processes
- Define auto-correlation function (ACF) and partial auto-correlation function (PACF) and use R to plot the sample ACF and sample PACF of a time series
- Explain the concepts of differencing and smoothing via moving averages to remove/highlight trends and seasonal components in a time series
- Define the zero-mean autoregressive process of order one or AR(1) and use R to obtain samples from this type of process
- Perform maximum likelihood estimation for the full and conditional likelihood in an AR(1)
- Perform Bayesian inference for the AR(1) under the conditional likelihood and the reference prior

**Week 2: The AR(p) process**

This module extends the concepts learned in Week 1 about the AR(1) process to the general case of the AR(p). Maximum likelihood estimation and Bayesian posterior inference in the AR(p) are discussed.

Learning Objectives
- Define the autoregressive process of order p or AR(p) and use R to obtain samples from such process
- Define ARIMA (autoregressive moving average) models (honors)
- Perform posterior inference for the AR(p) under the conditional likelihood and the reference prior
- Perform a full data analysis in R using an AR(p) including likelihood estimation and Bayesian inference, model order selection, and forecasting
- Explain the relationship between the AR characteristic polynomial, the ACF, the forecast function and the spectral density in the case of an AR(p)

**Week 3: Normal dynamic linear models, Part I**

Normal Dynamic Linear Models (NDLMs) are defined and illustrated in this module using several examples. Model building based on the forecast function via the superposition principle is explained. Methods for Bayesian filtering, smoothing and forecasting for NDLMs in the case of known observational variances and known system covariance matrices are discussed and illustrated.

Learning Objectives
- Use R for analysis and forecasting of time series using NDLM (case of known observational and system variances)
- Derive the equations to obtain posterior inference and forecasting in the NDLM with known observational and system variances, including the filtering, smoothing and forecasting equations
- Apply the NDLM superposition principle and explain the role of the forecast function
- Define trend and regression normal DLMs
- Explain the general normal dynamic linear model (NDLM) representation

**Week 4: Normal dynamic linear models, Part II**

Learning Objectives
- Use R for analysis and forecasting of time series using the NDLM (cases of known or unknown observational variance and unknown system variance specified using a discount factor)
- Derive the equations to obtain posterior inference and forecasting in the NDLM with unknown observational variance and system variance specified via discount factors
- Define seasonal NDLMs
- Apply the NDLM superposition principle and explain the role of the forecast function

**Week 5: Final Project**

In this final project you will use normal dynamic linear models to analyze a time series dataset downloaded from Google trend.

Learning Objectives
- Use R for analysis and forecasting of time series using NDLM (case of known observational and system variances)
- Use R for analysis and forecasting of time series using the NDLM (cases of known or unknown observational variance and unknown system variance specified using a discount factor)


## Course 5: Capstone Project

This is the capstone project for UC Santa Cruz's Bayesian Statistics Specialization. It is an opportunity for you to demonstrate a wide range of skills and knowledge in Bayesian statistics and to apply what you know to real-world data. You will review essential concepts in Bayesian statistics with lecture videos and quizzes, and you will perform a complex data analysis and compose a report on your methods and results.

**Week 1**



**Week 2**



**Week 3**



**Week 4**




# Certificate

![Certificate]()