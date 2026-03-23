data {
  int<lower=1> n_obs;
  vector[n_obs] obs;
  array[n_obs] int age;
  array[n_obs] int year;
  array[n_obs] int fleet;

  int<lower=1> na;
  int<lower=1> ny;
  int<lower=1> min_age;
  int<lower=1> min_year;

  matrix[na, ny] M;
  matrix[na, ny] stockMeanWeight;
  matrix[na, ny] propMature;
  real surveyTime;
}

parameters {
  vector[na] logN1Y;
  vector[ny] eps_std;
  vector[ny] logFY;
  vector[na] logFA;
  real logSdCatch;
  vector[na] logQ;
  real logSdSurvey;
  real logSdR;
  real tthetaR;
  real logMuR;
}

transformed parameters {
  matrix[na, ny] F;
  matrix[na, ny] logN;
  vector[ny] eps;
  real theta = 2 * inv_logit(tthetaR) - 1;
  real sdR = exp(logSdR);

  //Non-centered innovations
  eps = eps_std * sdR;

  //Fishing Mortality
  for (y in 1:ny) {
    for (a in 1:na) {
      F[a, y] = exp(logFA[a] + logFY[y]);
    }
  }

  //Population Dynamics with MA(1) Recruitment
  // Initial Year (Age-structure in year 1)
  for (a in 1:na) logN[a, 1] = logN1Y[a];

  for (y in 2:ny) {
    // Recruitment (Age 1) follows: logMuR + eps[y] + theta * eps[y-1]
    // RTMB loop used eps[i+1] + theta*eps[i]
    logN[1, y] = logMuR + eps[y] + theta * eps[y - 1];

    for (a in 2:na) {
      logN[a, y] = logN[a - 1, y - 1] - F[a - 1, y - 1] - M[a - 1, y - 1];
    }
  }
}

model {
  // Priors
  eps_std ~ std_normal(); // This is equivalent to eps ~ dnorm(0, sdR)
  logSdCatch ~ normal(0, 1);
  logSdSurvey ~ normal(0, 1);
  logSdR ~ normal(0, 1);
  // tthetaR ~ normal(0, 2);
  // logMuR ~ normal(10, 5);
  // logFY ~ normal(0, 2);
  // logFA ~ normal(0, 2);
  // logQ ~ normal(-5, 3);
  // logN1Y ~ normal(10, 5);

  // Likelihood
  for (i in 1:n_obs) {
    int a_idx = age[i] - min_age + 1;
    int y_idx = year[i] - min_year + 1;
    real logPred;
    real sd;

    if (fleet[i] == 1) { // Catch
      logPred = log(F[a_idx, y_idx]) - log(F[a_idx, y_idx] + M[a_idx, y_idx]) +
                log1m_exp(-F[a_idx, y_idx] - M[a_idx, y_idx]) + logN[a_idx, y_idx];
      sd = exp(logSdCatch);
    } else { // Survey
      logPred = logQ[a_idx] - (F[a_idx, y_idx] + M[a_idx, y_idx]) * surveyTime + logN[a_idx, y_idx];
      sd = exp(logSdSurvey);
    }

    log(obs[i]) ~ normal(logPred, sd);
  }
}

generated quantities {
  vector[ny] ssb;
  for (y in 1:ny) {
    ssb[y] = sum(rows_dot_product(exp(logN[, y]), (stockMeanWeight[, y] .* propMature[, y])));
  }
}
