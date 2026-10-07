# Conceptual figure

out_file <- "figures/sparsity_concept.png"
FIG_W <- 26.0
FIG_H <- 19.3

## Text sizes were tuned at 23in wide and the vertical gaps at 16.3in tall;
## scaling both keeps the layout proportional when the page size changes.
TEXT_BOOST <- 1.11      # raise for bigger type; the page grows to suit
HS <- FIG_W / 23.0 * TEXT_BOOST
VS <- 16.3 / FIG_H

## ---- greys ----------------------------------------------------------------
col_cell <- "#E0E0E0"
col_gry  <- "#565656"
col_lgry <- "#BDBDBD"
col_hdr  <- "#F0F0F0"
col_hdr2 <- "#DADADA"

CEX_PANEL <- 2.10 * HS
CEX_BODY  <- 1.72 * HS
CEX_SUB   <- 1.42 * HS
CEX_EQ    <- 1.62 * HS
CEX_STEP  <- 1.34 * HS
CEX_HOLD  <- 1.40 * HS
CEX_NODE  <- 1.25 * HS
CEX_LAB   <- 1.30 * HS
CEX_COST  <- 1.50 * HS
CEX_TERM  <- 1.45 * HS

## ---- geometry helpers -----------------------------------------------------
ipx   <- function() par("pin")[1] / diff(par("usr")[1:2])
ipy   <- function() par("pin")[2] / diff(par("usr")[3:4])
ry_of <- function(rx) rx * ipx() / ipy()

## ---- text helpers ---------------------------------------------------------
wrap_text <- function(txt, maxw, cex, font = 1) {
  words <- strsplit(txt, " ")[[1]]
  lines <- character(0); cur <- ""
  for (w in words) {
    test <- if (nchar(cur)) paste(cur, w) else w
    if (strwidth(test, cex = cex, font = font) > maxw && nchar(cur)) {
      lines <- c(lines, cur); cur <- w
    } else cur <- test
  }
  c(lines, cur)
}

lh_of <- function(cex) strheight("Ag", cex = cex) * 1.50

para <- function(x, y, txt, maxw, cex = CEX_BODY, col = "black", font = 1) {
  ls <- wrap_text(txt, maxw, cex, font)
  lh <- lh_of(cex)
  for (i in seq_along(ls))
    text(x, y - (i - 1) * lh, ls[i], adj = c(0, 1), cex = cex, col = col, font = font)
  invisible(y - length(ls) * lh)
}

## ---- drawing helpers ------------------------------------------------------
## The graph and the Hessian share a row inside a column, so both are sized
## from COLW rather than fixed, and the node spacing is whatever is left over.
NN      <- 5                       # nodes in the graph / rows in H
NODE_RX <- NULL                    # set once COLW is known
node <- function(x, y, lab) {
  ry <- ry_of(NODE_RX)
  th <- seq(0, 2 * pi, length.out = 72)
  polygon(x + NODE_RX * cos(th), y + ry * sin(th), col = "white",
          border = "black", lwd = 1.3)
  text(x, y, lab, cex = CEX_NODE)
}

arc <- function(x1, x2, y, depth, lwd = 1) {
  s <- seq(0, 1, length.out = 90)
  lines(x1 + (x2 - x1) * s, y - depth * sin(pi * s), col = "black", lwd = lwd)
}

spy <- function(x0, y0, w, n, nzfun) {
  cw <- w / n
  ch <- ry_of(cw)
  for (i in 1:n) for (j in 1:n)
    if (nzfun(i, j))
      rect(x0 + (j - 1) * cw, y0 - i * ch, x0 + j * cw, y0 - (i - 1) * ch,
           col = col_cell, border = "black", lwd = 0.8)
  rect(x0, y0 - n * ch, x0 + n * cw, y0, border = "black", lwd = 1.2)
  invisible(y0 - n * ch)
}

sep_rule <- function(y, x1 = 1.0, x2 = 99.0)
  lines(c(x1, x2), c(y, y), col = col_lgry, lwd = 1.1)

## Anything drawn inside a case column must fit it; collect violations rather
## than discovering them by eye in the rendered png.
.overflow <- character(0)
fits <- function(what, txt, cex, avail) {
  w <- strwidth(txt, cex = cex)
  if (w > avail)
    .overflow <<- c(.overflow, sprintf("%s: %.2f > %.2f units", what, w, avail))
  invisible(w)
}

## A step is either list(expression, grey) or a plain string shown as a grey
## note.  Heights are measured before anything is drawn so every column can put
## its graph, H and cost line on one shared baseline.
EQ_H <- 3.40 * VS
TERM_H <- 2.95 * VS
step_h <- function(steps, colw) {
  h <- 0
  for (st in steps)
    h <- h + if (is.character(st))
      length(wrap_text(st, colw - 3.2, CEX_STEP)) * lh_of(CEX_STEP) + 0.9 * VS
    else EQ_H
  h
}

## height of the trailing "the two likelihood terms" block, if a column has one
term_h <- function(term) {
  if (is.null(term)) 0 else 1.8 * VS + lh_of(CEX_STEP) +
    length(term[[2]]) * (lh_of(CEX_STEP) + 0.2 * VS + TERM_H)
}

## ===========================================================================
png(out_file, width = FIG_W, height = FIG_H, units = "in", res = 300)
op <- par(mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0), xpd = NA, family = "sans")
plot.new()
plot.window(xlim = c(0, 100), ylim = c(0, 100), xaxs = "i", yaxs = "i")

LM   <- 1.5
EW   <- 20.5                     # (e) is a full-height column on the right
ESEP <- 2.6                      # gap between the left region and (e)
EX   <- 100 - LM - EW            # left edge of (e)
RW   <- EX - ESEP - LM           # width available to (a)-(d)
AX   <- LM                       # (a) and the case columns start here

GUT  <- 1.3
GRP  <- 2.8
COLW <- (RW - 2 * GUT - GRP) / 4
NODE_RX <- COLW * 0.052
CX   <- c(LM,
          LM + COLW + GUT,
          LM + 2 * COLW + GUT + GRP,
          LM + 3 * COLW + 2 * GUT + GRP)

ytop <- 99.0

## ---------------------------------------------------------------------------
## (a)  where the cost comes from -- a wide, shallow band across the left
## ---------------------------------------------------------------------------
ATW <- 40.0                      # the prose column inside (a)
text(AX, ytop, "(a)  Where the cost comes from", adj = c(0, 1),
     cex = CEX_PANEL, font = 2)
y <- ytop - lh_of(CEX_PANEL) - 0.8

yl <- para(AX, y, paste(
  "The Laplace approximation integrates the random effects u out of the joint negative",
  "log-likelihood f: minimise f over u, then add half the log determinant of its Hessian H.",
  "Every gradient step therefore factorizes H, at a cost determined entirely by which entries",
  "of H are nonzero. That is a property of which variable was declared random."),
  ATW, CEX_BODY)

EQA <- expression(-log ~ L(theta) %~~% f(theta, hat(u)) +
                  textstyle(frac(1, 2)) ~ log ~ det ~ H(theta, hat(u)))
text(AX + ATW + 2.5, y - 1.4, EQA, adj = c(0, 1), cex = 1.62 * HS)
fits("panel a eq1", EQA, 1.62 * HS, 24.0)
text(AX + ATW + 2.5, y - 6.4, "up to an additive constant", adj = c(0, 1),
     cex = 1.30 * HS, col = col_gry, font = 3)

EQB <- expression(H[jk] == frac(partialdiff^2 * f,
                                partialdiff * u[j] * partialdiff * u[k]))
text(AX + ATW + 27.0, y - 1.4, EQB, adj = c(0, 1), cex = 1.62 * HS)
fits("panel a eq2", EQB, 1.62 * HS, RW - ATW - 27.0)

y <- min(yl, y - 10.0) - 1.6
RULE <- expression(paste(u[j], " and ", u[k], " are joined ") %<->%
                   paste(" some term of f contains both"))
rect(AX, y - 4.6, AX + RW, y, col = col_hdr, border = col_lgry)
text(AX + RW / 2, y - 2.3, RULE, adj = c(0.5, 0.5), cex = 1.30 * HS)
ya <- y - 5.6

## ---------------------------------------------------------------------------
## the right-hand column: (b) on top, then (e) below a dividing rule
## ---------------------------------------------------------------------------
lines(rep(EX - ESEP / 2, 2), c(ytop + 0.6, 1.0), col = col_lgry, lwd = 1.1)

ECEX <- CEX_PANEL * 0.88
PLOTH <- 17.0                    # height reserved for the cost plot

text(EX, ytop, "(b)  Why sparsity is faster", adj = c(0, 1), cex = ECEX, font = 2)
fits("panel b title", "(b)  Why sparsity is faster", ECEX, EW)
plot_top <- ytop - lh_of(ECEX) - 0.8
plot_bot <- plot_top - PLOTH

CAP <- paste("TMB factorizes H with a sparse Cholesky: it only visits entries that",
             "can be nonzero, so zeros cost nothing.")
ye <- para(EX, plot_bot - 1.6, CAP, EW, 1.46 * HS, col = col_gry)

## the small rule that separates the two panels in this column
ye <- ye - 2.4
lines(c(EX, EX + EW), c(ye, ye), col = col_lgry, lwd = 1.1)
ye <- ye - 2.8

text(EX, ye, "(e)  The plus group", adj = c(0, 1), cex = ECEX, font = 2)
fits("panel e title", "(e)  The plus group", ECEX, EW)
ye <- ye - lh_of(ECEX) - 0.6

ye <- para(EX, ye, paste(
  "Here the recruitment deviations are the random effects and everything else is",
  "deterministic survival. Even so, the oldest age keeps its own survivors. Write",
  "P for the plus group, C for the cohort entering it, and S for survival:"),
  EW, CEX_BODY, col = col_gry)
ye <- ye - 1.6

for (st in list(
  list(expression(P[t] == C[t] * S + P[t-1] * S), FALSE),
  list(expression(P[t-1] == C[t-1] * S + P[t-2] * S), TRUE),
  "the same line, one year back",
  list(expression(P[t] == C[t] * S + C[t-1] * S^2 + P[t-2] * S^2), FALSE),
  "a plus group term is still on the right, so the step runs again, and again",
  list(expression(P[t] == C[t] * S + C[t-1] * S^2 + C[t-2] * S^3 + ldots), FALSE))) {
  if (is.character(st)) {
    ye <- para(EX + 1.0, ye, st, EW - 2.0, CEX_STEP, col = col_gry) - 0.9 * VS
  } else {
    text(EX + 0.8, ye, st[[1]], adj = c(0, 1), cex = CEX_EQ,
         col = if (isTRUE(st[[2]])) col_gry else "black")
    fits("panel e eq", st[[1]], CEX_EQ, EW - 1.6)
    ye <- ye - EQ_H
  }
}
ye <- ye - 1.0

ye <- para(EX, ye, paste(
  "every C is a different cohort, and every cohort carries one recruitment",
  "deviation, so the plus group holds all of them back to year one and H is dense"),
  EW, CEX_HOLD, font = 2)
ye <- ye - 2.4

text(EX, ye, "two ways to break it:", adj = c(0, 1),
     cex = CEX_LAB, col = col_gry, font = 3)
ye <- ye - lh_of(CEX_LAB) - 1.0

fixes <- list(
  list("remove the plus group",
       "the oldest age becomes a true age, so nothing feeds into itself"),
  list("auxiliary variables",
       "give ages 2+ their own process error (SAM, WHAM); the state space grows"))
for (fx in fixes) {
  text(EX + 0.6, ye, fx[[1]], adj = c(0, 1), cex = CEX_TERM, font = 2)
  fits("plusgroup label", fx[[1]], CEX_TERM, EW - 1.2)
  ye <- ye - lh_of(CEX_TERM) - 0.4
  ye <- para(EX + 1.4, ye, fx[[2]], EW - 2.0, CEX_STEP, col = col_gry) - 1.4
}

## band 1 ends where panel (a) does; (e) runs the full height beside everything
band1_bot <- ya - 1.2
sep_rule(band1_bot, x1 = LM, x2 = EX - ESEP)

## ===========================================================================
## one case column
## ===========================================================================
## `steps` is a list of list(expr, grey); `notes` is grey text shown where a
## sparse column's substitution would have continued.  `ygraph` is passed in
## so every column puts its graph, H and cost line on the same baseline.
draw_case <- function(x0, ytop, ygraph, dense, headline, lab, steps,
                      holds, nodelab, cost, term = NULL) {

  rect(x0, ytop - 2.9 * VS, x0 + COLW, ytop,
       col = if (dense) col_hdr2 else col_hdr, border = "black", lwd = 1.0)
  text(x0 + 1.2, ytop - 1.45 * VS, headline, adj = c(0, 0.5), cex = CEX_SUB, font = 2)
  fits(paste("header", headline), headline, CEX_SUB, COLW - 2.4)
  y <- ytop - 4.3 * VS

  text(x0 + 1.0, y, lab, adj = c(0, 1), cex = CEX_LAB, col = col_gry, font = 3)
  fits(paste("lab", headline), lab, CEX_LAB, COLW - 2.0)
  y <- y - lh_of(CEX_LAB) - 0.8 * VS

  for (st in steps) {
    if (is.character(st)) {
      y <- para(x0 + 1.8, y, st, COLW - 3.2, CEX_STEP, col = col_gry) - 0.9 * VS
    } else {
      text(x0 + 1.3, y, st[[1]], adj = c(0, 1), cex = CEX_EQ,
           col = if (isTRUE(st[[2]])) col_gry else "black")
      fits(paste("eq", headline), st[[1]], CEX_EQ, COLW - 2.3)
      y <- y - EQ_H
    }
  }

  ## the likelihood term this buys, where a column has the room for it
  if (!is.null(term)) {
    y <- y - 0.9 * VS
    text(x0 + 1.8, y, term[[1]], adj = c(0, 1), cex = CEX_STEP, col = col_gry,
         font = 3)
    fits(paste("termlab", headline), term[[1]], CEX_STEP, COLW - 3.6)
    y <- y - lh_of(CEX_STEP) - 0.9 * VS
    for (tr in term[[2]]) {
      text(x0 + 1.7, y, tr[[1]], adj = c(0, 1), cex = CEX_STEP, col = col_gry)
      y <- y - lh_of(CEX_STEP) - 0.2 * VS
      text(x0 + 2.4, y, tr[[2]], adj = c(0, 1), cex = CEX_TERM)
      fits(paste("term", headline), tr[[2]], CEX_TERM, COLW - 3.2)
      y <- y - TERM_H
    }
  }

  if (y < ygraph + 5.4 * VS)
    .overflow <<- c(.overflow,
                    sprintf("VERTICAL %s: text reaches %.1f, baseline is %.1f",
                            headline, y, ygraph + 5.4 * VS))

  ## the punchline sits on a shared baseline in all four columns
  text(x0 + 1.0, ygraph + 5.4 * VS, holds, adj = c(0, 1), cex = CEX_HOLD, font = 2)
  fits(paste("holds", headline), holds, CEX_HOLD, COLW - 2.0)

  nn  <- NN
  mw  <- COLW * 0.27
  mx0 <- x0 + COLW - mw - 0.8
  gx0 <- x0 + 1.0 + NODE_RX
  gsp <- ((mx0 - 0.8) - (x0 + 1.0) - 2 * NODE_RX) / (nn - 1)

  text(x0 + 1.0, ygraph, "dependency graph", adj = c(0, 1), cex = CEX_LAB,
       col = col_gry)
  text(mx0, ygraph, "Hessian H", adj = c(0, 1), cex = CEX_LAB, col = col_gry)
  ytop2 <- ygraph - lh_of(CEX_LAB) - 0.6 * VS

  mtop <- ytop2
  grow <- mtop - ry_of(mw) / 2
  dmax <- 4.5 * VS
  if (dense) {
    for (i in 1:(nn - 1)) for (j in (i + 1):nn)
      arc(gx0 + (i - 1) * gsp, gx0 + (j - 1) * gsp, grow - 1.6 * VS,
          depth = 1.2 + (dmax - 1.2) * ((j - i) / (nn - 1))^0.8, lwd = 0.9)
  } else {
    for (i in 1:(nn - 1))
      lines(c(gx0 + (i - 1) * gsp, gx0 + i * gsp), c(grow, grow),
            col = "black", lwd = 2.1)
  }
  for (i in 1:nn) node(gx0 + (i - 1) * gsp, grow, nodelab(i))

  nzf  <- if (dense) function(i, j) TRUE else function(i, j) abs(i - j) <= 1
  mbot <- spy(mx0, mtop, mw, nn, nzf)

  yy <- min(grow - 1.6 * VS - dmax, grow - ry_of(NODE_RX) - 0.6 * VS, mbot) - 1.4 * VS
  text(x0 + 1.0, yy - lh_of(CEX_COST), cost, adj = c(0, 0), cex = CEX_COST,
       font = 2)
  fits(paste("cost", headline), cost, CEX_COST, COLW - 2.0)
  invisible(yy - lh_of(CEX_COST))
}


# c and d -----------------------------------------------------------------


y <- band1_bot - 2.0
ytitle <- y

text(CX[1], y, expression(paste(bold("(c)  AR(1):"), "    ",
                                bolditalic(x[t]) == bolditalic(phi) *
                                bolditalic(x[t-1]) + bolditalic(epsilon[t]))),
     adj = c(0, 1), cex = CEX_PANEL)
text(CX[3], y, expression(paste(bold("(d)  MA(1):"), "    ",
                                bolditalic(x[t]) == bolditalic(epsilon[t]) +
                                bolditalic(theta) * bolditalic(epsilon[t-1]))),
     adj = c(0, 1), cex = CEX_PANEL)
y <- y - lh_of(CEX_PANEL) - 0.6

ARW <- CX[2] + COLW - CX[1]
yb <- para(CX[1], y, paste(
  "Each likelihood term has to be written using only the random effects. That rewrite finishes in",
  "one step if the states are random, but never finishes if the innovations are."),
  ARW, 1.58 * HS, col = col_gry)
yc <- para(CX[3], y, paste(
  "The same rewrite with the roles swapped. Now it finishes in one step if the innovations are",
  "random, and never finishes if the states are."),
  ARW, 1.58 * HS, col = col_gry)
y <- min(yb, yc) - 1.8
ytop_cases <- y

## ---- the four substitution chains -----------------------------------------
NOTE_STOP <- paste("Both quantities on the right are already random effects, so there is",
                   "nothing more to replace and the rewrite is done.")

st_ar_sparse <- list(
  list(expression(epsilon[t] == x[t] - phi * x[t-1]), FALSE),
  NOTE_STOP)

st_ar_dense <- list(
  list(expression(x[t] == phi * x[t-1] + epsilon[t]), FALSE),
  list(expression(x[t-1] == phi * x[t-2] + epsilon[t-1]), TRUE),
  "the same line, one year back",
  list(expression(x[t] == phi * (phi * x[t-2] + epsilon[t-1]) + epsilon[t]), FALSE),
  list(expression(x[t] == phi^2 * x[t-2] + phi * epsilon[t-1] + epsilon[t]), FALSE),
  paste("The right still contains a state, so it has to be replaced too, bringing in",
        "the year before, and so on back to year one."),
  list(expression(x[t] == epsilon[t] + phi * epsilon[t-1] +
                          phi^2 * epsilon[t-2] + ldots), FALSE))

st_ma_sparse <- list(
  list(expression(x[t] == epsilon[t] + theta * epsilon[t-1]), FALSE),
  NOTE_STOP)

st_ma_dense <- list(
  list(expression(epsilon[t] == x[t] - theta * epsilon[t-1]), FALSE),
  list(expression(epsilon[t-1] == x[t-1] - theta * epsilon[t-2]), TRUE),
  "the same line, one year back",
  list(expression(epsilon[t] == x[t] - theta * (x[t-1] - theta * epsilon[t-2])), FALSE),
  list(expression(epsilon[t] == x[t] - theta * x[t-1] + theta^2 * epsilon[t-2]), FALSE),
  paste("The right still contains an innovation, so it has to be replaced too, bringing",
        "in the year before, and so on back to year one."),
  list(expression(epsilon[t] == x[t] - theta * x[t-1] +
                                theta^2 * x[t-2] - ldots), FALSE))

HEAD2 <- "year t terms:"

tm_ar_sparse <- list(HEAD2, list(
  list("process:",     expression(-log ~ N(x[t] ~ "|" ~ phi * x[t-1], sigma))),
  list("observation:", expression(-log ~ N(Y[t] ~ "|" ~ x[t], tau)))))

tm_ar_dense <- list(HEAD2, list(
  list("process:",     expression(-log ~ N(epsilon[t] ~ "|" ~ 0, sigma))),
  list("observation:", expression(-log ~ N(Y[t] ~ "|" ~ Sigma[k <= t] ~
                                           phi^{t-k} * epsilon[k], tau)))))

tm_ma_sparse <- list(HEAD2, list(
  list("process:",     expression(-log ~ N(epsilon[t] ~ "|" ~ 0, sigma))),
  list("observation:", expression(-log ~ N(Y[t] ~ "|" ~ epsilon[t] +
                                           theta * epsilon[t-1], tau)))))

tm_ma_dense <- list(HEAD2, list(
  list("process:",     expression(-log ~ N(Sigma[k <= t] ~ (-theta)^{t-k} * x[k] ~
                                           "|" ~ 0, sigma))),
  list("observation:", expression(-log ~ N(Y[t] ~ "|" ~ x[t], tau)))))

## one shared baseline, set by the tallest chain
YGRAPH <- ytop_cases - 4.3 * VS - lh_of(CEX_LAB) - 0.8 * VS -
  max(mapply(function(st, tm) step_h(st, COLW) + term_h(tm),
             list(st_ar_sparse, st_ar_dense, st_ma_sparse, st_ma_dense),
             list(tm_ar_sparse, tm_ar_dense, tm_ma_sparse, tm_ma_dense))) - 6.0 * VS

k1 <- draw_case(
  CX[1], ytop_cases, YGRAPH, dense = FALSE,
  "random effect: states", "rewrite \u03b5 using x alone:", st_ar_sparse,
  expression(paste("process term: ", group("{", list(x[t-1], x[t]), "}"))),
  function(i) bquote(x[.(i)]),
  expression(paste(bold("SPARSE"), "   \u00b7   ", 2 * n, " operations")),
  term = tm_ar_sparse)

k2 <- draw_case(
  CX[2], ytop_cases, YGRAPH, dense = TRUE,
  "random effect: innovations", "rewrite x using \u03b5 alone:", st_ar_dense,
  expression(paste("observation term: ", group("{", list(epsilon[1], ldots, epsilon[t]), "}"))),
  function(i) bquote(epsilon[.(i)]),
  expression(paste(bold("DENSE"), "   \u00b7   ", n^3 / 3, " operations")),
  term = tm_ar_dense)

k3 <- draw_case(
  CX[3], ytop_cases, YGRAPH, dense = FALSE,
  "random effect: innovations", "rewrite x using \u03b5 alone:", st_ma_sparse,
  expression(paste("observation term: ", group("{", list(epsilon[t-1], epsilon[t]), "}"))),
  function(i) bquote(epsilon[.(i)]),
  expression(paste(bold("SPARSE"), "   \u00b7   ", 2 * n, " operations")),
  term = tm_ma_sparse)

k4 <- draw_case(
  CX[4], ytop_cases, YGRAPH, dense = TRUE,
  "random effect: states", "rewrite \u03b5 using x alone:", st_ma_dense,
  expression(paste("process term: ", group("{", list(x[1], ldots, x[t]), "}"))),
  function(i) bquote(x[.(i)]),
  expression(paste(bold("DENSE"), "   \u00b7   ", n^3 / 3, " operations")),
  term = tm_ma_dense)

## vertical divider between the AR pair and the MA pair
lines(rep((CX[2] + COLW + CX[3]) / 2, 2),
      c(ytitle + 1.0, min(k1, k2, k3, k4) - 0.5), col = col_lgry, lwd = 1.1)

cat(sprintf("COLW=%.1f  a_bot=%.1f  e_bot=%.1f  band1_bot=%.1f  cases_bot=%.1f\n",
            COLW, ya, ye, band1_bot, min(k1, k2, k3, k4)))
if (length(.overflow)) {
  cat("OVERFLOW:\n"); cat(paste0("  ", .overflow, collapse = "\n"), "\n")
} else cat("no overflow\n")
stopifnot(min(k1, k2, k3, k4) > 0.5, ye > 0.5)

## ---- the cost plot, drawn last so it owns its coordinate system ---------
par(fig = c(EX / 100, (EX + EW) / 100, plot_bot / 100, plot_top / 100),
    new = TRUE, xpd = FALSE, mar = c(3.2, 4.6, 0.6, 0.6), mgp = c(1.6, 0.6, 0))
nn <- seq(20, 600, length.out = 400)
plot(nn, nn^3 / 3, type = "n", log = "xy", xlab = "", ylab = "",
     axes = FALSE, ylim = c(30, 1e8), xaxs = "i")
lines(nn, nn^3 / 3, col = "black", lwd = 3.4)
lines(nn, 2 * nn,   col = "black", lwd = 3.4, lty = 2)
box(col = "black", lwd = 1.2)
axis(1, at = c(20, 50, 100, 200, 500), cex.axis = 1.35 * HS, col = "black",
     col.axis = "black", tcl = -0.35, lwd = 1.1)
axis(2, at = 10^seq(2, 8, by = 2), labels = parse(text = paste0("10^", seq(2, 8, by = 2))),
     cex.axis = 1.35 * HS, col = "black", col.axis = "black", las = 1,
     tcl = -0.35, lwd = 1.1)
mtext("years in the series (n)", side = 1, line = 2.0, cex = 1.35 * HS)
mtext("operations", side = 2, line = 3.3, cex = 1.35 * HS)
text(70, 25 * 70^3 / 3, "dense",  cex = 1.5 * HS, font = 2, adj = c(0.5, 0))
text(330, 3 * 2 * 330,   "sparse", cex = 1.5 * HS, font = 2, adj = c(0.5, 0))

par(op)
dev.off()
cat("wrote", out_file, "\n")
