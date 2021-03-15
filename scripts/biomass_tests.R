## ---- libraries
library(knitr)
library(tidyverse)
library(simstudy)
library(glmmTMB)
library(DHARMa)
library(emmeans)
library(plotly)
library(crosstalk)
library(INLA)
library(INLAutils)
library(brms)
library(broom.mixed)
library(rstan)
library(inlatools)
library(tidybayes)
library(mgcv)
library(gratia)
## ----end


## ---- dirs
if (!dir.exists('../data/')) dir.create('../data')
## ----end


## ---- synthetic-data-1
set.seed(123)

## Reef level
def.reef <- defData(varname='Reef', dist='normal', formula=0, variance=2.5, id='idReef')
def.reef <- def.reef %>% defData(varname='nSites', dist='nonrandom', formula=3)
dtReef <- genData(6, def.reef)

dtSite <- dtReef %>% genCluster("idReef", numIndsVar = "nSites", level1ID = "idSite")
def.site <- defData(varname='Site', dist='normal',formula=0, variance=0.5)
def.site <- def.site %>% defData(varname='nTransects', dist='nonrandom', formula=5)
dtSite <- def.site %>% addColumns(dtSite)

dtTransect <- dtSite %>% genCluster("idSite", numIndsVar = "nTransects", level1ID = "idTransect")
def.transect <- defData(varname="Transect", dist="normal", formula=0, variance=2.5)
#def.transect <- def.transect %>% defData(varname="nGroups", dist="noZeroPoisson", formula = 8)
dtTransect <- def.transect %>% addColumns(dtTransect)
dtTransect
## ----end

## ---- synthetic-data-2
## Create a different temporal trend per Species

theta = cbind(c(0.6,0.9,0.8,0.3,0.8,0.5,0.2),
              c(0.2,0.5,0.4,0.3,0.8,0.5,0.4),
              c(0.2,0.3,0.4,0.5,0.3,0.5,0.3),
              c(0.6,0.8,0.8,0.3,0.4,0.4,0.3),
              c(0.5,0.5,0.4,0.3,0.9,0.6,0.4),
              c(0.8,0.9,0.7,0.5,0.5,0.3,0.1),
              c(0.9,0.5,0.4,0.3,0.1,0.1,0.1),
              c(0.7,0.6,0.7,0.5,0.4,0.4,0.2)
              )
knots=c(0.25,0.5,0.75)
## visualise these trends
viewSplines(theta=theta, knots=knots, degree=3)
## ----end

## ---- synthetic-data-3
dtYear <- addPeriods(dtTransect, nPeriods=10, idvars='idTransect') %>%
    genSpline(newvar='Species1',predictor='period', theta=theta[,1], knots=knots, noise.var=0.0, newrange='0;1') %>%
    genSpline(newvar='Species2',predictor='period', theta=theta[,2], knots=knots, noise.var=0.0, newrange='0;1') %>% 
    genSpline(newvar='Species3',predictor='period', theta=theta[,3], knots=knots, noise.var=0.0, newrange='0;1') %>% 
    genSpline(newvar='Species4',predictor='period', theta=theta[,4], knots=knots, noise.var=0.0, newrange='0;1') %>% 
    genSpline(newvar='Species5',predictor='period', theta=theta[,5], knots=knots, noise.var=0.0, newrange='0;1') %>% 
    genSpline(newvar='Species6',predictor='period', theta=theta[,6], knots=knots, noise.var=0.0, newrange='0;1') %>% 
    genSpline(newvar='Species7',predictor='period', theta=theta[,7], knots=knots, noise.var=0.0, newrange='0;1') %>% 
    genSpline(newvar='Species8',predictor='period', theta=theta[,8], knots=knots, noise.var=0.0, newrange='0;1') %>%
    pivot_longer(cols=c(Species1, Species2, Species3, Species4, Species5, Species5, Species6, Species7, Species8),
                 names_to = 'Species') %>%
    mutate(Species=as.numeric(gsub('Species([0-9])','\\1',Species))) %>%
    data.table:::as.data.table()
dtYear
dtYear %>%
    ggplot() +
    geom_line(aes(y=value, x=period, color=as.factor(Species)))
## ----end

## ---- synthetic-data-4
dtGroup <- dtYear %>%
    group_by(idReef,Reef,idSite,Site,idTransect,Transect, period) %>%
    sample_n(size=rpois(1,6), replace=TRUE, weight=1:8/sum(1:8)) %>%
    mutate(idGroup=1:n()) %>%
    ungroup

dtGroup %>% group_by(Species) %>% count()
## ----end

## ---- synthetic-data-5
def.data = defDataAdd(varname='Actual_Abundance',dist='poisson', formula='3+(value*2)+Reef+Site+Transect+Species/2')
def.data <- def.data %>% defData(varname='zi',dist='nonrandom', formula=0)
def.data <- def.data %>% defData(varname='Abundance',formula='zi | .40 + Actual_Abundance | .60', dist='mixture')
def.data <- def.data %>% defDataAdd(varname='Length', dist='gamma', formula="50/Species", variance=1)


dat <- addColumns(def.data, dtGroup %>% data.table::as.data.table()) %>%
    mutate(TRANSECT=factor(idTransect),
           SITE=factor(idSite),
           REEF=factor(idReef),
           YEAR=period,
           fYEAR=factor(period),
           SPECIES=factor(Species),
           GROUP=factor(idGroup),
           ABUNDANCE=Abundance,
           LENGTH=Length) %>%
    dplyr::select(REEF, SITE, TRANSECT, YEAR, fYEAR, SPECIES, GROUP, ABUNDANCE, LENGTH)

dat %>%
    ggplot() +
    geom_line(aes(y=ABUNDANCE, x=YEAR, color=SPECIES, group=interaction(TRANSECT,GROUP))) +
    facet_wrap(~SITE)
dat %>%
    ggplot() +
    geom_line(aes(y=LENGTH, x=YEAR, color=SPECIES, group=interaction(TRANSECT,GROUP))) +
    scale_y_sqrt() + 
    facet_wrap(~SITE)
## ----end



## ---- interactive_fig1
dat.sum=dat %>%
    group_by(SITE, TRANSECT, YEAR, SPECIES) %>%
    summarise(ABUNDANCE=sum(ABUNDANCE,na.rm=TRUE),
              LENGTH=mean(LENGTH, na.rm=TRUE)) %>% # aggregate over groups of fish
    ungroup

sites = unique(dat.sum$SITE)
data.shared <- SharedData$new(dat.sum)
g1 = data.shared %>%
    ggplot() +
    geom_line(aes(y=ABUNDANCE, x=YEAR, color=SPECIES, group=TRANSECT))
g2 = data.shared %>%
    ggplot() +
    geom_line(aes(y=LENGTH, x=YEAR, color=SPECIES, group=TRANSECT))
 
widgets <- bscols(filter_select("Sites", "Select sites", data.shared, ~SITE)
                  )
                                        #widgets <- bscols(filter_checkbox("Sites", "Sites", dd, ~idSite, columns=3))
gp <- subplot(ggplotly(g1) %>% style(showlegend=FALSE), ggplotly(g2),
        shareX = TRUE, titleY=TRUE,
        margin=0.1,
        which_layout=2
        )
bscols(widths=c(2,8), widgets,gp
       )
## ----end

## ---- interactive_fig2
g1 = data.shared %>%
    ggplot() +
    geom_line(aes(y=ABUNDANCE, x=YEAR, color=SPECIES, group=TRANSECT)) +
    facet_wrap(~SITE)
ggplotly(g1)
## ----end


## Data processing #########################################################

## ---- data-processing
lw.lookup <- tribble(
    ~SPECIES, ~a, ~b,
    1,0.029,2.3,
    2,0.013,2.5,
    3,0.009,3.1,
    4,0.025,3.0,
    5,0.027,2.9,
    6,0.037,2.9,
    7,0.032,2.7,
    8,0.030,2.1
) %>% mutate(SPECIES=factor(SPECIES))
lw.lookup

dat = dat %>%
    left_join(lw.lookup) %>%
    mutate(Biomass=a * ABUNDANCE * LENGTH^b) %>%
    group_by(REEF, SITE, TRANSECT, YEAR, fYEAR) %>%
    summarise(ABUNDANCE=sum(ABUNDANCE),
              BIOMASS=sum(Biomass)/1000) %>%
    ungroup %>%
    dplyr::select(REEF, SITE, TRANSECT, YEAR, fYEAR, ABUNDANCE, BIOMASS)
dat

dat %>%
    ggplot() +
    geom_line(aes(y=ABUNDANCE, x=fYEAR, group=TRANSECT)) +
    facet_wrap(~REEF)

dat %>%
    ggplot() +
    geom_line(aes(y=BIOMASS, x=fYEAR, group=TRANSECT)) +
    facet_wrap(~REEF) +
    scale_y_log10()

save(dat, file='../data/dat.RData')
## ----end

## Modelling options ###############################################################

## Abundance -----------------------------------------------------------------------

## ZINB (glmmTMB) ..................................................................

## ---- fitModel_zinb
mod <- glmmTMB(ABUNDANCE ~ fYEAR + (1|REEF)+(1|SITE)+(1|TRANSECT),
               zi = ~fYEAR,
               data=dat, family=nbinom2)
## ----end

## ---- fitModel_zinb_diagnostics
resids <- DHARMa::simulateResiduals(mod, plot=TRUE)
testZeroInflation(resids)
## ----end

## ---- fitModel_zinb_summary
summary(mod)
## ----end

## ---- fitModel_zinb_fig
emmeans(mod, ~fYEAR, type='response') %>%
    as.data.frame %>%
    ggplot(aes(y=response, x=as.numeric(fYEAR))) +
    geom_ribbon(aes(ymin=lower.CL, ymax=upper.CL), alpha=0.3) +
    geom_line()
## ----end

## ZINB (INLA) ..................................................................

## ---- fitModel_zinb_INLA
newdata = with(dat,
               data.frame(fYEAR=levels(fYEAR),
                          REEF=NA, SITE=NA, TRANSECT=NA, ABUNDANCE=NA))
dat.pred = dat %>% bind_rows(newdata)
zinb_inla.mod <- inla(ABUNDANCE ~ fYEAR + f(REEF, model='iid')+f(SITE, model='iid')+
                          f(TRANSECT, model='iid'),
                      data=dat.pred, family='zeroinflatednbinomial1',
                     # control.family=list(link='log'),
                      control.predictor=list(link=1, compute=TRUE),
                      control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE, config=TRUE)
                      )
save(zinb_inla.mod, file='../data/zinb_inla.mod.RData')
## ----end

## ---- fitModel_zinb_INLA_posteriors
load(file='../data/zinb_inla.mod.RData')
preds <- lapply(1:nrow(dat), function(x) inla.rmarginal(250, zinb_inla.mod$marginals.fitted.values[[x]]))
#preds <- lapply(1:nrow(dat), function(x) inla.rmarginal(250, zinb_inla.mod$marginals.linear.predictor[[x]]))
preds = do.call('cbind', preds)
#preds = exp(preds)
#preds[1:10,1:10]
resids <- createDHARMa(simulatedResponse = t(preds),
                             observedResponse = dat$ABUNDANCE,
                             fittedPredictedResponse = apply(preds, 2, median),
                             integerResponse = TRUE)
plot(resids)

autoplot(zinb_inla.mod)
ggplot_inla_residuals(zinb_inla.mod, dat$ABUNDANCE)
ggplot_inla_residuals2(zinb_inla.mod, dat$ABUNDANCE)

fast_distribution_check(zinb_inla.mod) %>% plot()
# ratio of observed to expected remains close to 1 (100%) - in particular the envelope includes the reference of 100%
dispersion_check(zinb_inla.mod) %>% plot
## If there is no under or over dispersion P(D|data < D|model) ~ 0.5).
## If it is close to 0 = underdispersed.
## If it is close to 1 = over dispersed
## Dashed line indicates dispersion of original data
## P(D|data > D|model) in title
## Above seems underdispersed 
## ----end

## ---- fitModel_zinb_INLA_summary
summary(zinb_inla.mod)
## ----end

## ---- fitModel_zinb_INLA_fig
newdata <- cbind(newdata,
                 zinb_inla.mod$summary.fitted.values[(nrow(dat)+1):nrow(dat.pred),]) %>%
    rename(lower=`0.025quant`, upper=`0.975quant`, median=`0.5quant`)
newdata %>%
    ggplot(aes(y=median, x=as.numeric(fYEAR))) +
    geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.3) +
    geom_line()
## ----end

## ZINB (brms) ..................................................................

## ---- fitModel_zinb_brm
mod.brm <- brm(bf(ABUNDANCE ~ fYEAR + (1|REEF)+(1|SITE)+(1|TRANSECT),
                  zi = ~fYEAR,
                  family='zero_inflated_negbinomial'),
               ## prior=c(
               ##     prior(normal(1,0.1), class='Intercept'),
               ##     prior(normal(0,0.1), class='b'),
               ##     prior(cauchy(0,1), class='sd'),
               ##     prior(gamma(0.01,0.01), class='shape'),
               ##     prior(logistic(0,1), class="Intercept", dpar='zi'),
               ##     prior(normal(0,0.1), class='b', dpar='zi')
               ##     ),
               data=dat, 
               iter=3000,
               warmup=1000,
               thin=5,
               chains=3, cores=3,
               control = list(adapt_delta=0.99)
               )
save(mod.brm, file='../data/mod.brm.RData')
## ----end

## ---- fitModel_zinb_brm_MCMCdiagnostics
load(file='../data/mod.brm.RData')
stan_trace(mod.brm$fit)
stan_ac(mod.brm$fit)
stan_rhat(mod.brm$fit)
stan_ess(mod.brm$fit)
## ----end

## ---- fitModel_zinb_brm_diagnostics
preds <- posterior_predict(mod.brm,  nsamples=250,  summary=FALSE)
mod.resids <- createDHARMa(simulatedResponse = t(preds),
                           observedResponse = dat$ABUNDANCE,
                           fittedPredictedResponse = apply(preds, 2, median),
                           integerResponse = TRUE)
plot(mod.resids)
testZeroInflation(mod.resids)
## ----end

## ---- fitModel_zinb_brm_summary
summary(mod.brm)
## ----end

## ---- fitModel_zinb_brm_fig
load(file='../data/mod.brm.RData')
emmeans(mod.brm, ~fYEAR, type='response') %>%
    as.data.frame %>%
    ggplot(aes(y=prob, x=as.numeric(fYEAR))) +
    geom_ribbon(aes(ymin=lower.HPD, ymax=upper.HPD), alpha=0.3) +
    geom_line()
## ----end

## Biomass -------------------------------------------------------------------------

## Tweedie (INLA) ..................................................................

## ---- fitModel_tweedie_INLA
newdata = with(dat,
               data.frame(fYEAR=levels(fYEAR),
                          REEF=NA, SITE=NA, TRANSECT=NA, BIOMASS=NA))
dat.pred = dat %>% bind_rows(newdata)

tweedie_inla.mod <- inla(BIOMASS ~ fYEAR + f(REEF, model='iid')+f(SITE, model='iid')+
                          f(TRANSECT, model='iid'),
                         data=dat.pred,
                         family='tweedie',
                         control.family=list(link='log'),
                         control.predictor=list(link=1, compute=TRUE),
                         control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE, config=TRUE)
                         )
## ----end

## Delta method (INLA) .............................................................

## ---- fitModel_detla_INLA
draws=1000
newdata = with(dat,
               data.frame(fYEAR=levels(fYEAR),
                          REEF=NA, SITE=NA, TRANSECT=NA, BIOMASS=NA))
dat.pred = dat %>% bind_rows(newdata)
## binomial model
delta_inla.mod.b <- inla(BIOMASS ~ fYEAR + f(REEF, model='iid')+f(SITE, model='iid')+
                          f(TRANSECT, model='iid'),
                         data=dat.pred %>% mutate(BIOMASS=as.numeric(BIOMASS>0)),
                         family='binomial',
                         control.family=list(link='logit'),
                         #control.predictor=list(link=1, compute=TRUE),
                         control.predictor=list(link=NA, compute=TRUE),
                         control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE, config=TRUE)
                         )

save(delta_inla.mod.b, file='../data/delta_inla.mod.b.RData')
b = newdata %>%
    mutate(ID=1:n()+nrow(dat)) %>%
    ## slice(-1:-5) %>%
    dplyr::select(-BIOMASS) %>%
    group_by(fYEAR) %>%
    nest() %>%
    mutate(Fit.b=map(.x=data, .f=function(.x) inla.rmarginal(draws,delta_inla.mod.b$marginals.fitted.values[[.x$ID]]))) %>%
    #mutate(Fit.b=map(.x=data, .f=function(.x) inla.rmarginal(draws,mod.inla.b$marginals.linear.predictor[[.x$ID]])),
    #       Fit.b=plogis(Fit.b)) %>%
    dplyr::select(-data) %>%
    unnest(cols=c(Fit.b)) %>%
    mutate(Iter=1:n(),
           Fit.b=plogis(Fit.b)) %>%
    ungroup()

## Gamma
dat.g = dat %>% filter(BIOMASS!=0) %>% droplevels
dat.pred = dat.g %>% bind_rows(newdata)
delta_inla.mod.g <- inla(BIOMASS ~ fYEAR + f(REEF, model='iid')+f(SITE, model='iid')+
                          f(TRANSECT, model='iid'),
                         data=dat.pred,
                         family='gamma',
                         control.family=list(link='log'),
                         control.predictor=list(link=NA, compute=TRUE),
                         control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE, config=TRUE)
                         )
g = newdata %>%
    mutate(ID=1:n()+nrow(dat.g)) %>%
    ## slice(-1:-5) %>%
    dplyr::select(-BIOMASS) %>%
    group_by(fYEAR) %>%
    nest() %>%
    #mutate(Fit.g=map(.x=data, .f=function(.x) inla.rmarginal(draws,mod.inla.g$marginals.fitted.values[[.x$ID]]))) %>%
    mutate(Fit.g=map(.x=data, .f=function(.x) inla.rmarginal(draws, delta_inla.mod.g$marginals.linear.predictor[[.x$ID]]))) %>%
    dplyr::select(-data) %>%
    unnest(cols=c(Fit.g)) %>%
    mutate(Iter=1:n(),
           Fit.g=exp(Fit.g)) %>%
    #mutate(Iter=1:n()) %>%
    ungroup()

## Combine together
newdata <- full_join(newdata,
                     b %>% full_join(g) %>%
                     mutate(Fit=exp(log(Fit.b) + log(Fit.g))) %>%
                     group_by(fYEAR) %>%
                     median_qi(Fit)) %>%
    dplyr::rename(median=Fit,lower=.lower, upper=.upper) %>%
    dplyr::select(-.point,-.interval, -.width)
## ----end

## ---- fitModel_delta_INLA_fig
newdata %>%
    ggplot(aes(y=median, x=as.numeric(fYEAR))) +
    geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.3) +
    geom_line()
## ----end

## Tweedie (mgcv) ..................................................................

## ---- fitModel_tweedie_mgcv
tweedie_mgcv.mod <- gam(BIOMASS ~ fYEAR + s(REEF, bs='re')+s(SITE, bs='re')+
                          s(TRANSECT, bs='re'), method='REML',
                        data=dat,
                        family='tw')
## ----end

## ---- fitModel_tweedie_mgcv_diagnostics
gratia::appraise(tweedie_mgcv.mod)
## ----end

## ---- fitModel_tweedie_mgcv_summary
summary(tweedie_mgcv.mod)
## ----end

## ---- fitModel_tweedie_mgcv_fig
emmeans(tweedie_mgcv.mod, ~fYEAR, type='response') %>%
    as.data.frame %>%
    ggplot(aes(y=response, x=as.numeric(fYEAR))) +
    geom_ribbon(aes(ymin=lower.CL, ymax=upper.CL), alpha=0.3) +
    geom_line()
## ----end

## Tweedie (glmmTMB) ...............................................................

## ---- fitModel_tweedie_glmmTMB
tweedie_glmmTMB.mod <- glmmTMB(BIOMASS ~ fYEAR + (1|REEF)+(1|SITE)+(1|TRANSECT),
               data=dat, family=tweedie)
save(tweedie_glmmTMB.mod, file='../data/tweedie_glmmTMB.mod.RData')
## ----end

## ---- fitModel_tweedie_glmmTMB_diagnostics
resids <- DHARMa::simulateResiduals(tweedie_glmmTMB.mod, plot=TRUE)
testZeroInflation(resids)
## ----end

## ---- fitModel_tweedie_glmmTMB_summary
summary(tweedie_glmmTMB.mod)
## ----end

## ---- fitModel_tweedie_glmmTMB_fig
emmeans(tweedie_glmmTMB.mod, ~fYEAR, type='response') %>%
    as.data.frame %>%
    ggplot(aes(y=response, x=as.numeric(fYEAR))) +
    geom_ribbon(aes(ymin=lower.CL, ymax=upper.CL), alpha=0.3) +
    geom_line()
## ----end

## giGamma (glmmTMB) ..............................................................

## ---- fitModel_zigamma_glmmTMB
zigamma_glmmTMB.mod <- glmmTMB(BIOMASS ~ fYEAR + (1|REEF)+(1|SITE)+(1|TRANSECT),
               zi = ~fYEAR,
               data=dat, family=ziGamma)
save(zigamma_glmmTMB.mod, file='../data/tweedie_glmmTMB.mod.RData')
## ----end

## ---- fitModel_zigamma_glmmTMB_diagnostics
#resids <- DHARMa::simulateResiduals(zigamma_glmmTMB.mod, plot=TRUE)
#testZeroInflation(resids)
## ----end

## ---- fitModel_zigamma_glmmTMB_summary
summary(zigamma_glmmTMB.mod)
## ----end

## ---- fitModel_zigamma_glmmTMB_fig
emmeans(zigamma_glmmTMB.mod, ~fYEAR, type='response') %>%
    as.data.frame %>%
    ggplot(aes(y=response, x=as.numeric(fYEAR))) +
    geom_ribbon(aes(ymin=lower.CL, ymax=upper.CL), alpha=0.3) +
    geom_line()
## ----end
