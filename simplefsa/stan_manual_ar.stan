data {
  int<lower=1> n_obs;
  vector[n_obs] obs;
  array[n_obs] int age;
  array[n_obs] int year;
  array[n_obs] int fleet;

  int<lower=1> na; // Number of ages
  int<lower=1> ny; // Number of years
  int<lower=1> min_age;
  int<lower=1> min_year;

  matrix[na, ny] M;
  matrix[na, ny] stockMeanWeight;
  matrix[na, ny] propMature;
  real surveyTime;
}

parameters {
  vector[na] logN1Y;
  vector[ny-1] logN1A;
  vector[ny] logFY;
  vector[na] logFA;
  real logSdCatch;
  vector[na] logQ; // Assuming logQ is age-specific for fleet 2
  real logSdSurvey;
  real logSdR;
  real tphiR;
  real logMuR;
  vector[ny - 1] z;
}

transformed parameters {
  matrix[na, ny] F;
  matrix[na, ny] logN;
  vector[ny - 1] rec_dev;
  real phi = 2 * inv_logit(tphiR) - 1;
  real sdR = exp(logSdR);

  // Setup Fishing Mortality
  for (y in 1:ny) {
    for (a in 1:na) {
      F[a, y] = exp(logFA[a] + logFY[y]);
    }
  }

  // Recruitment AR1 Process (Non-centered)
  rec_dev[1] = z[1] * (sdR / sqrt(1 - phi^2));
  for (y in 2:(ny - 1)) {
    rec_dev[y] = phi * rec_dev[y - 1] + z[y] * sdR;
  }

  // Population Dynamics Matrix
  // Initial Year
  for (a in 1:na) logN[a, 1] = logN1Y[a];

  for (y in 2:ny) {
    // Recruitment (Age 1)
    logN[1, y] = logMuR + rec_dev[y - 1];
    for (a in 2:na) {
      logN[a, y] = logN[a - 1, y - 1] - F[a - 1, y - 1] - M[a - 1, y - 1];
    }
  }
}

model {
  // Priors (Adjust based on your stock's biology)
  z ~ std_normal();
  logSdCatch ~ normal(0, 1);
  logSdSurvey ~ normal(0, 1);
  logSdR ~ normal(0, 1);
  // tphiR ~ normal(0, 2);
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

    if (fleet[i] == 1) { // Catch Fleet
      logPred = log(F[a_idx, y_idx]) - log(F[a_idx, y_idx] + M[a_idx, y_idx]) +
                log1m_exp(-F[a_idx, y_idx] - M[a_idx, y_idx]) + logN[a_idx, y_idx];
      sd = exp(logSdCatch);
    } else { // Survey Fleet
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
