<!DOCTYPE html>

<html xmlns="http://www.w3.org/1999/xhtml">

<head>

<meta charset="utf-8" />
<meta http-equiv="Content-Type" content="text/html; charset=utf-8" />
<meta name="generator" content="pandoc" />

<meta name="viewport" content="width=device-width, initial-scale=1">

<meta name="author" content="Mark Puttick" />

<meta name="date" content="2017-10-11" />





</head>

<body>




<h1 class="title toc-ignore">motmot</h1>
<h4 class="author"><em>Mark Puttick</em></h4>
<h4 class="date"><em>2017-10-11</em></h4>



<p>Models Of Trait Macroevolution On Trees (MOTMOT) is an R package that allows for testing of models of trait evolution <span class="citation">(Thomas and Freckleton 2012)</span>.</p>
<ul>
<li><a href="#models-of-trait-evolution">Tree transformation</a> models estimated using Maximum likelihood: <a href="#brownian-motion">Brownian motion</a>, <a href="pagels-lambda">Pagel’s Lambda</a>, <a href="#delta">Delta</a>, <a href="#kappa">Kappa</a>, <a href="#ornstein-uhlenbeck">Ornstein-Uhlenbeck (OU)</a>, <a href="#acdc">Acceleration-Deaceleration (ACDC)</a>, and <a href="#estimate-lambda-alongside-models">estimating lambda alongside other models</a></li>
<li><a href="#rate-heterogeneous-models-of-evolution">Rate heterogeneous models of evolution</a>. Fit models in which the rate of evolution differs in clades selected <a href="#rate-heterogeneity-selected-a-priori"><em>a priori</em></a> <span class="citation">(Thomas, Freckleton, and Székely 2006, <span class="citation">O’Meara et al. (2006)</span>)</span>, and models with <a href="#rate-heterogeneity-with-no-a-priori-information">no <em>a-priori</em> shift locations</a> <span class="citation">(Thomas and Freckleton 2012)</span></li>
<li><a href="#time-slice-model">TimeSlice</a> fit models in which all rates change at a specific time(s) by tested all times or those selected by the user</li>
<li><a href="#nested-models-of-evolution">Nested Shift mode</a> Fit models models in which the ancestral BM rate switches to a ‘nested’ rate within a monophyletic clade in the phylogeny</li>
<li><a href="#bayesian-estimation-of-tree-transformation-models">Bayesian estimation</a> of tree transformation models</li>
</ul>
<div id="introduction" class="section level1">
<h1>Introduction</h1>
<p>First we will install motmot.2.0 from github</p>
<div class="sourceCode"><pre class="sourceCode r"><code class="sourceCode r"><span class="cf">if</span>(<span class="op">!</span><span class="kw">any</span>(<span class="st">&quot;devtools&quot;</span> <span class="op">%in%</span><span class="st"> </span><span class="kw">rownames</span>(<span class="kw">installed.packages</span>()))) <span class="kw">install.packages</span>(<span class="st">&quot;devtools&quot;</span>)
<span class="kw">library</span>(devtools)
<span class="kw">install_github</span>(<span class="st">&quot;ghthomas/motmot&quot;</span>, <span class="dt">ref=</span><span class="st">&quot;motmot.2.0&quot;</span>)</code></pre></div>
<pre><code>## Downloading GitHub repo ghthomas/motmot@motmot.2.0
## from URL https://api.github.com/repos/ghthomas/motmot/zipball/motmot.2.0</code></pre>
<pre><code>## Installing motmot.2.0</code></pre>
<pre><code>## '/Library/Frameworks/R.framework/Resources/bin/R' --no-site-file  \
##   --no-environ --no-save --no-restore --quiet CMD INSTALL  \
##   '/private/var/folders/lb/grtql71j7d1d8wyln9l7x8prddz3_y/T/Rtmp3tKOGS/devtools2ed85aa4229/ghthomas-motmot-7ad80f7'  \
##   --library='/Users/mp1728/Library/R/3.4/library' --install-tests</code></pre>
<pre><code>## </code></pre>
<p>For these examples we will use anolis data available from motmot. A time-calibrated phylogeny of anolis species (“anolis.tree”), and various trait and biogeographical trait data (“anolis.data”)</p>
<div class="sourceCode"><pre class="sourceCode r"><code class="sourceCode r"><span class="kw">library</span>(motmot.<span class="fl">2.0</span>, <span class="dt">quietly=</span>T)
<span class="kw">data</span>(anolis.tree)
<span class="kw">data</span>(anolis.data)

<span class="kw">names</span>(anolis.data)</code></pre></div>
<pre><code>## [1] &quot;Species&quot;      &quot;Island_type&quot;  &quot;ecomorph&quot;     &quot;geo_ecomorph&quot;
## [5] &quot;Female_SVL&quot;   &quot;Male_SVL&quot;</code></pre>
<div class="sourceCode"><pre class="sourceCode r"><code class="sourceCode r"><span class="kw">attach</span>(anolis.data)
anolis.tree</code></pre></div>
<pre><code>## 
## Phylogenetic tree with 165 tips and 164 internal nodes.
## 
## Tip labels:
##  A_occultus, A_darlingt, A_monticol, A_bahoruco, A_dolichoc, A_henderso, ...
## Node labels:
##  2, 2, 2, 2, 2, 2, ...
## 
## Rooted; includes branch lengths.</code></pre>
<p>For the first part of the tutorial we will use the continuous trait data: male snout-ventral length ‘Male_SVL’. We will construct a matrix of just these data, and check if we have missing data</p>
<div class="sourceCode"><pre class="sourceCode r"><code class="sourceCode r">male.length &lt;-<span class="st"> </span><span class="kw">matrix</span>(Male_SVL, <span class="dt">dimnames=</span><span class="kw">list</span>(<span class="kw">rownames</span>(anolis.data)))
<span class="kw">any</span>(<span class="kw">is.na</span>(male.length[,<span class="dv">1</span>]))</code></pre></div>
<pre><code>## [1] TRUE</code></pre>
<p>We do. So we will remove these data from the male.length data, and log the trait data</p>
<div class="sourceCode"><pre class="sourceCode r"><code class="sourceCode r">complete.male.length &lt;-<span class="st"> </span><span class="kw">complete.cases</span>(male.length)
missing.species &lt;-<span class="st"> </span><span class="kw">rownames</span>(male.length)[<span class="op">!</span>complete.male.length]
male.length &lt;-<span class="st"> </span><span class="kw">as.matrix</span>(male.length[complete.male.length, ])
male.length &lt;-<span class="st"> </span><span class="kw">log</span>(male.length)</code></pre></div>
<p>Finally, we will ‘prune’ the species from the tree using ‘drop.tip’ from APE. Do our species from the data and tree now match?</p>
<div class="sourceCode"><pre class="sourceCode r"><code class="sourceCode r">phy &lt;-<span class="st"> </span><span class="kw">drop.tip</span>(anolis.tree, missing.species)
<span class="kw">name.check</span>(phy, male.length)</code></pre></div>
<pre><code>## [1] &quot;OK&quot;</code></pre>
</div>
<p>They do. We can now plot our tree and data using the “traitData.plot” function</p>
<div class="sourceCode"><pre class="sourceCode r"><code class="sourceCode r"><span class="kw">traitData.plot</span>(<span class="dt">y=</span>male.length, phy)</code></pre></div>
<p>They do. We can now plot our tree and data using the “traitData.plot” function</p>
<img atl=" " src="https://github.com/ghthomas/motmot/blob/motmot.2.0/vignettes/figures/plot1-1.png"/>
<p class="caption">traitData showing the realtive male snout ventral length at the tips</p>
<div id="models-of-trait-evolution" class="section level1">
<h1>Models of trait evolution</h1>
<p>We can now test various models of evolution using our trait data.</p>
<div id="brownian-motion" class="section level2">
<h2>Brownian motion</h2>
<p>To start we will fit a simple Brownian motion model to the data</p>
<div class="sourceCode"><pre class="sourceCode r"><code class="sourceCode r">bm.ml &lt;-<span class="st"> </span><span class="kw">transformPhylo.ML</span>(phy, <span class="dt">y=</span>male.length, <span class="dt">model=</span><span class="st">&quot;bm&quot;</span>)
bm.ml</code></pre></div>
<pre><code>## $brownianVariance
##             [,1]
## [1,] 0.001917681
## 
## $logLikelihood
## [1] -17.35662
## 
## $root.state
## [1] 4.120935
## 
## $AIC
## [1] 38.71324
## 
## $AICc
## [1] 38.78872</code></pre>
</div>
<div id="pagels-lambda" class="section level2">
<h2>Pagel’s lambda</h2>
<p>We can also fit models to test Pagel’s lambda</p>
<div class="sourceCode"><pre class="sourceCode r"><code class="sourceCode r">lambda.ml &lt;-<span class="st"> </span><span class="kw">transformPhylo.ML</span>(phy, <span class="dt">y=</span>male.length, <span class="dt">model=</span><span class="st">&quot;lambda&quot;</span>)
lambda.ml</code></pre></div>
<pre><code>## $MaximumLikelihood
## [1] -13.95609
## 
## $Lambda
##       MLLambda  LowerCI   UpperCI
## [1,] 0.9685786 0.879023 0.9970324
## 
## $brownianVariance
##             [,1]
## [1,] 0.001675265
## 
## $root.state
## [1] 4.122186
## 
## $AIC
## [1] 33.91218
## 
## $AICc
## [1] 34.06408</code></pre>
<p>We can see from these results the CI of lambda is &lt; 1, but there is still a large phylogenetic signal is these data</p>
<p>A new feature in motmot allows for plotting of the likelihood profile for the branch-transformation parameter, in this case Pagel’s lambda</p>
<div class="sourceCode"><pre class="sourceCode r"><code class="sourceCode r">lambda.ml &lt;-<span class="st"> </span><span class="kw">transformPhylo.ML</span>(phy, <span class="dt">y=</span>male.length, <span class="dt">model=</span><span class="st">&quot;lambda&quot;</span>, <span class="dt">profilePlot=</span>T)</code></pre></div>
<div class="figure">
<img src="https://github.com/ghthomas/motmot/blob/motmot.2.0/vignettes/figures/plot2-1.png" alt=" " />
<p class="caption">profile plot of ML estimation for Pagel’s lambda</p>
</div>
<p>We can now compare the fit of the BM and Lambda models. Lambda has higher likelihood, but it also has more parameters. We can test whether this is a significant improvement. First we will use the chi-squared distribution. The models differ in one degree of freedom: BM has 2 parameters (brownian variance, root state) and lambda has those two plus the value of lambda. We can use the R function pchisq to obtain a p value, and see that lambda is not a superior fit to these male length data</p>
<div class="sourceCode"><pre class="sourceCode r"><code class="sourceCode r">p.value &lt;-<span class="st"> </span><span class="dv">1</span> <span class="op">-</span><span class="st"> </span><span class="kw">pchisq</span>(lambda.ml<span class="op">$</span>MaximumLikelihood <span class="op">-</span><span class="st"> </span>bm.ml<span class="op">$</span>logLikelihood, <span class="dv">1</span>)
p.value</code></pre></div>
<pre><code>## [1] 0.06517542</code></pre>
<p>However, there is a large Akaike Information Criterion (AICc) difference between the two models: BM has a higher AICc compared to Lambda. The differce (4.724636) is &gt;4 which is tradtionally seen as indication of a superior fit (Burnham and Anderson 2003).</p>
<div class="sourceCode"><pre class="sourceCode r"><code class="sourceCode r">bm.ml<span class="op">$</span>AICc<span class="op">-</span><span class="st"> </span>lambda.ml<span class="op">$</span>AICc</code></pre></div>
<pre><code>## [1] 4.724636</code></pre>
<p>The parameters, brownian variance, root state, Maximum likelihoods, AIC, and AICc can be obtained for a number of models in motmot.</p>
</div>
<div id="delta" class="section level2">
<h2>Delta</h2>
<p>Delta indicates a slow or increase in the rate of trait evolution through time; a value of 1 is equivalent to Brownian motion, &lt; 1 indicates a slow-down, and &gt; 1 is difficult to interpret (greater change near the present). Here we find a MLE of 1.27 but the CI spans &lt; 1</p>
<div class="sourceCode"><pre class="sourceCode r"><code class="sourceCode r">delta.ml &lt;-<span class="st"> </span><span class="kw">transformPhylo.ML</span>(phy, <span class="dt">y=</span>male.length, <span class="dt">model=</span><span class="st">&quot;delta&quot;</span>, <span class="dt">profilePlot=</span>T)</code></pre></div>
<div class="figure">
<img src="https://github.com/ghthomas/motmot/blob/motmot.2.0/vignettes/figures/plot3-1.png" alt="" />
<p class="caption">profile plot to estimate delta</p>
</div>
<div class="sourceCode"><pre class="sourceCode r"><code class="sourceCode r">delta.ml</code></pre></div>
<pre><code>## $MaximumLikelihood
## [1] -17.05532
## 
## $Delta
##      MLDelta   LowerCI  UpperCI
## [1,] 1.27924 0.6092799 2.127126
## 
## $brownianVariance
##              [,1]
## [1,] 0.0004598011
## 
## $root.state
## [1] 4.127836
## 
## $AIC
## [1] 40.11064
## 
## $AICc
## [1] 40.26254</code></pre>
</div>
<div id="kappa" class="section level2">
<h2>Kappa</h2>
<p>Kappa is used as a measure of punctuated evolution and spans values of 0-1. 1 is equivalent to BM, and 0 indicates trait change occurs at events of speciation. Here the is some evidence of punctuated evolution, but the CI spans a value greater than one (as seen in the warning message)</p>
<div class="sourceCode"><pre class="sourceCode r"><code class="sourceCode r">kappa.ml &lt;-<span class="st"> </span><span class="kw">transformPhylo.ML</span>(phy, <span class="dt">y=</span>male.length, <span class="dt">model=</span><span class="st">&quot;kappa&quot;</span>, <span class="dt">profilePlot=</span>T)</code></pre></div>
<pre><code>## Warning in transformPhylo.ML(phy, y = male.length, model = &quot;kappa&quot;,
## profilePlot = T): Confidence limits fall outside the parameter bounds -
## consider changing lowerBound and/or upperBound</code></pre>
<div class="figure">
<img src="https://github.com/ghthomas/motmot/blob/motmot.2.0/vignettes/figures/plot4-1.png" alt="" />
<p class="caption">profile plot to estimate kappa</p>
</div>
</div>
<div id="ornstein-uhlenbeck" class="section level2">
<h2>Ornstein-Uhlenbeck</h2>
<p>The OU model allows for modelling of attraction to a optimum value (alpha), and psi fits a acceleration-deacceleration model to assess to the relative contributions of speciation and gradual evolution to a trait’s evolutionary rate</p>
<div class="sourceCode"><pre class="sourceCode r"><code class="sourceCode r">ou.ml &lt;-<span class="st"> </span><span class="kw">transformPhylo.ML</span>(phy, <span class="dt">y=</span>male.length, <span class="dt">model=</span><span class="st">&quot;OU&quot;</span>, <span class="dt">profilePlot=</span>T)</code></pre></div>
<pre><code>## [1] &quot;Warning - different start values produces different OU likelihoods&quot;</code></pre>
<pre><code>## Warning in transformPhylo.ML(phy, y = male.length, model = &quot;OU&quot;,
## profilePlot = T): Confidence limits fall outside parameter bounds -
## consider changing lowerBound and/or upperBound</code></pre>
<div class="figure">
<img src="https://github.com/ghthomas/motmot/blob/motmot.2.0/vignettes/figures/plot5-1.png" alt="profile plot to estimate alpha" />
<p class="caption">profile plot to estimate alpha</p>
</div>
<div class="sourceCode"><pre class="sourceCode r"><code class="sourceCode r">ou.ml</code></pre></div>
<pre><code>## $MaximumLikelihood
## [1] -16.84778
## 
## $Alpha
##          MLAlpha LowerCI    UpperCI
## [1,] 0.003690604      NA 0.01144679
## 
## $brownianVariance
##             [,1]
## [1,] 0.002206995
## 
## $root.state
## [1] 4.124566
## 
## $AIC
## [1] 39.69557
## 
## $AICc
## [1] 39.84747</code></pre>
</div>
<div id="acdc" class="section level2">
<h2>ACDC</h2>
<p>A new addition to MOTMOT is the ACDC model <span class="citation">(Harmon et al. 2010; Blomberg, Jr, and Ives 2003)</span>. This model allows for exponential changes in the rate of evolution in the history of a clade. If the upperBound value is &lt; 0, this is equivalent to the ‘Early Burst’ model fit in geiger</p>
<div class="sourceCode"><pre class="sourceCode r"><code class="sourceCode r">acdc.ml &lt;-<span class="st"> </span><span class="kw">transformPhylo.ML</span>(phy, <span class="dt">y=</span>male.length, <span class="dt">model=</span><span class="st">&quot;ACDC&quot;</span>, <span class="dt">profilePlot=</span>T)</code></pre></div>
<div class="figure">
<img src="https://github.com/ghthomas/motmot/blob/motmot.2.0/vignettes/figures/plot6-1.png" alt="" />
<p class="caption">profile plot to estimate the ACDC parameter</p>
</div>
<div class="sourceCode"><pre class="sourceCode r"><code class="sourceCode r">acdc.ml</code></pre></div>
<pre><code>## $MaximumLikelihood
## [1] -16.84778
## 
## $ACDC
##           MLacdc      LowerCI    UpperCI
## [1,] 0.007384034 -0.006756934 0.02288467
## 
## $brownianVariance
##             [,1]
## [1,] 0.001117947
## 
## $root.state
## [1] 4.124568
## 
## $AIC
## [1] 39.69557
## 
## $AICc
## [1] 39.84747</code></pre>
<p>There is little evidence here of exponential decreases or increases in the rate of trait evolution - the acdc exponential parameter is close to 0 (0.007). We can see this is not a significant improvement on BM</p>
<div class="sourceCode"><pre class="sourceCode r"><code class="sourceCode r">p.value.<span class="dv">2</span> &lt;-<span class="st"> </span><span class="dv">1</span> <span class="op">-</span><span class="st"> </span><span class="kw">pchisq</span>(acdc.ml<span class="op">$</span>MaximumLikelihood <span class="op">-</span><span class="st"> </span>bm.ml<span class="op">$</span>logLikelihood , <span class="dv">1</span>)
p.value.<span class="dv">2</span></code></pre></div>
<pre><code>## [1] 0.4756421</code></pre>
</div>
<div id="estimate-lambda-alongside-models" class="section level2">
<h2>Estimate lambda alongside models</h2>
<p>One way to deal with ‘noisy’ data is to estimate Pagel’s lambda alongside a parameter of interest. In motmot, lambda can be estimated alongside the delta, kappa, OU, psi, and ACDC models. Here we look at example using ACDC. The model is fit with same function. ‘transformPhyo.ML’, but with the argument ‘lambdaEst’ set to TRUE</p>
<div class="sourceCode"><pre class="sourceCode r"><code class="sourceCode r">acdc.ml.lambda &lt;-<span class="st"> </span><span class="kw">transformPhylo.ML</span>(<span class="dt">y=</span>male.length, phy, <span class="dt">model=</span><span class="st">&quot;ACDC&quot;</span>, <span class="dt">lambdaEst=</span>T)
<span class="co"># original ACDC model</span>
acdc.ml</code></pre></div>
<pre><code>## $MaximumLikelihood
## [1] -16.84778
## 
## $ACDC
##           MLacdc      LowerCI    UpperCI
## [1,] 0.007384034 -0.006756934 0.02288467
## 
## $brownianVariance
##             [,1]
## [1,] 0.001117947
## 
## $root.state
## [1] 4.124568
## 
## $AIC
## [1] 39.69557
## 
## $AICc
## [1] 39.84747</code></pre>
<div class="sourceCode"><pre class="sourceCode r"><code class="sourceCode r"><span class="co"># ACDC model plus lambda</span>
acdc.ml.lambda</code></pre></div>
<pre><code>## $MaximumLikelihood
## [1] -13.89171
## 
## $ACDC
##            MLacdc    LowerCI    UpperCI
## [1,] -0.003593448 -0.0176518 0.01111201
## 
## $brownianVariance
##             [,1]
## [1,] 0.002101792
## 
## $root.state
## [1] 4.121082
## 
## $lambda
## [1] 0.9583652
## 
## $AIC
## [1] 37.78342
## 
## $AICc
## [1] 38.16804</code></pre>
<p>We can see lambda is &gt; 1, and this has affected the parameter estimation (slightly). However, the improvement in the model fit is not significant compared to the ACDC model without lambda, or the null BM model</p>
<div class="sourceCode"><pre class="sourceCode r"><code class="sourceCode r"><span class="co"># p value of the ACDC and ACDC+lambda models. No significant improvement</span>
<span class="dv">1</span> <span class="op">-</span><span class="st"> </span><span class="kw">pchisq</span>(acdc.ml.lambda<span class="op">$</span>MaximumLikelihood <span class="op">-</span><span class="st"> </span>acdc.ml<span class="op">$</span>MaximumLikelihood , <span class="dt">df=</span><span class="dv">1</span>)</code></pre></div>
<pre><code>## [1] 0.08555553</code></pre>
<div class="sourceCode"><pre class="sourceCode r"><code class="sourceCode r"><span class="co"># p value of the BM and ACDC+lambda model comparison. No significant improvement</span>
<span class="dv">1</span> <span class="op">-</span><span class="st"> </span><span class="kw">pchisq</span>(acdc.ml.lambda<span class="op">$</span>MaximumLikelihood <span class="op">-</span><span class="st"> </span>bm.ml<span class="op">$</span>logLikelihood, <span class="dt">df=</span><span class="dv">2</span>)</code></pre></div>
<pre><code>## [1] 0.1768496</code></pre>
</div>
</div>
<div id="rate-heterogeneous-models-of-evolution" class="section level1">
<h1>Rate heterogeneous models of evolution</h1>
<div id="rate-heterogeneity-selected-a-priori" class="section level2">
<h2>rate heterogeneity selected <em>a priori</em></h2>
<p>MOTMOT can test models of evolution in which pre-defined clades can vary in the rate of evolution. Here we fit a model in which the nodes descending from nodes 183 and 240 have a seperate rate of evolution. We can visualise these nodes on the phylogeny</p>
<div class="sourceCode"><pre class="sourceCode r"><code class="sourceCode r"><span class="kw">plot</span>(phy, <span class="dt">show.tip.label=</span>F, <span class="dt">no.margin=</span>T, <span class="dt">edge.col=</span><span class="st">&quot;grey20&quot;</span>)
<span class="kw">nodelabels</span>(<span class="kw">c</span>(<span class="dv">183</span>, <span class="dv">240</span>), <span class="kw">c</span>(<span class="dv">183</span>, <span class="dv">240</span>), <span class="dt">bg=</span><span class="st">&quot;black&quot;</span>, <span class="dt">col=</span><span class="st">&quot;white&quot;</span>)</code></pre></div>
<div class="figure">
<img src="https://github.com/ghthomas/motmot/blob/motmot.2.0/vignettes/figures/plot7-1.png" alt="" />
<p class="caption">lineages with different rates of evolution</p>
</div>
<p>We then fit the motmot model, again using the function transformPhylo.ML. We use the argument “model=clade”. This fits the non-censored model of O’Meara et al. (2006).</p>
<div class="sourceCode"><pre class="sourceCode r"><code class="sourceCode r">cladeRate.ml &lt;-<span class="st"> </span><span class="kw">transformPhylo.ML</span>(male.length, <span class="dt">phy=</span>phy, <span class="dt">model=</span><span class="st">&quot;clade&quot;</span>, <span class="dt">nodeIDs=</span><span class="kw">c</span>(<span class="dv">183</span>, <span class="dv">240</span>))
cladeRate.ml</code></pre></div>
<pre><code>## $Rates
##      node    MLRate   LowerCI UpperCI
## [1,]  183 0.7853043 0.3223851 2.47607
## [2,]  240 0.8119501 0.3927771 2.01348
## 
## $MaximumLikelihood
## [1] -17.1441
## 
## $brownianVariance
##             [,1]
## [1,] 0.001969926
## 
## $root.state
## [1] 4.123463
## 
## $AIC
## [1] 42.28819
## 
## $AICc
## [1] 42.54297</code></pre>
<p>These results indicate that the two clades tend to have a lower rate of evolution compared to the background rate. However, the CIs indicate these decreases may not be robust</p>
</div>
<div id="rate-heterogeneity-with-no-a-priori-information" class="section level2">
<h2>rate heterogeneity with no <em>a priori</em> information</h2>
<p>We can also fit rate heterogeneous models without specifying where we expect shifts on the tree. We can use the arguments “model=”tm1“” and “model=”tm2“”; these models fit ‘traitMedusa’ models in which all nodes are tested for rate increases or decreases. It is possible to exclude small nodes using the argument ‘minCladeSize’. As well as allowing clade differences in rate, the “tm2” also allows for branch-based increases or decreases in rate.</p>
<p>We will fit the “tm2” model that allows for clade- and branch-specific changes in rate. This uses the familiar function ‘transformPhylo.ML’. As the ‘traitMedusa’ is computationally-intensive (fitting models to all nodes, optimising increasing numbers of rates) we will fit the ‘tm2’ to a subset of the phylogeny and data using the APE function ‘extract.clade’</p>
<div class="sourceCode"><pre class="sourceCode r"><code class="sourceCode r"><span class="kw">plot</span>(phy, <span class="dt">show.tip.label=</span>F, <span class="dt">no.margin=</span>T, <span class="dt">edge.col=</span><span class="st">&quot;grey20&quot;</span>)
<span class="kw">nodelabels</span>(<span class="dv">182</span>, <span class="dv">182</span>, <span class="dt">bg=</span><span class="st">&quot;black&quot;</span>, <span class="dt">col=</span><span class="st">&quot;white&quot;</span>)</code></pre></div>
<div class="figure">
<img src="https://github.com/ghthomas/motmot/blob/motmot.2.0/vignettes/figures/plot8-1.png" alt="" />
<p class="caption">the subset of the tree</p>
</div>
<div class="sourceCode"><pre class="sourceCode r"><code class="sourceCode r">phy.clade &lt;-<span class="st"> </span><span class="kw">extract.clade</span>(phy, <span class="dv">182</span>)
male.length.clade &lt;-<span class="st"> </span><span class="kw">as.matrix</span>(male.length[<span class="kw">match</span>(phy.clade<span class="op">$</span>tip.label, <span class="kw">rownames</span>(male.length)),])</code></pre></div>
<p>We can now fit the ‘tm2’ algorithm. The output shows the log-likelihood, AIC, AICc, rate type (branch of clade), for the best-fitting model at each stage. This starts with the BM model, and then one shift model, two shift model, etc.,</p>
<div class="sourceCode"><pre class="sourceCode r"><code class="sourceCode r"><span class="co"># not run</span>
<span class="co"># tm1.ml &lt;- transformPhylo.ML(y=male.length.clade, phy=phy.clade, model=&quot;tm1&quot;, minCladeSize=2, nSplits=3)</span>
<span class="co"># trait.medusa.tm1.summary &lt;- traitMedusaSummary(tm1.ml, cutoff=2, AICc=T)</span>
tm2.ml &lt;-<span class="st"> </span><span class="kw">transformPhylo.ML</span>(<span class="dt">y=</span>male.length.clade, <span class="dt">phy=</span>phy.clade, <span class="dt">model=</span><span class="st">&quot;tm2&quot;</span>, <span class="dt">minCladeSize=</span><span class="dv">2</span>, <span class="dt">nSplits=</span><span class="dv">3</span>)</code></pre></div>
<pre><code>## 
##  BM model
##         node shiftPos        lnL n.params      AIC     AICc
## BM         0        1 -0.2838554        2 4.567711 5.047711
## 
##  Shift 1
##         node shiftPos    lnL n.params       AIC       AICc   rate.1
## shift.1   45    clade 3.6472        3 -1.294399 -0.2943992 8.603397
## 
##  Shift 2
##         node shiftPos      lnL n.params       AIC      AICc   rate.1
## shift.2   39    clade 6.109191        5 -2.218382 0.5088912 7.278441
##           rate.2
## shift.2 0.137804
## 
##  Shift 3
##         node shiftPos      lnL n.params       AIC     AICc  rate.1
## shift.3   31    clade 8.628804        7 -3.257608 2.342392 6.68661
##            rate.2     rate.3
## shift.3 0.1261663 0.00541799</code></pre>
<p>We can now summarise the results of these data using ‘traitMedusaSummary’ and plotting the shifts on the phylogeny using ‘plotPhylo.motmot’. These results show an increase at node 45 that we can visualise on the phylogeny.</p>
<div class="sourceCode"><pre class="sourceCode r"><code class="sourceCode r">trait.medusa.tm2.summary &lt;-<span class="st"> </span><span class="kw">traitMedusaSummary</span>(tm2.ml, <span class="dt">cutoff=</span><span class="dv">2</span>, <span class="dt">AICc=</span>T)
trait.medusa.tm2.summary</code></pre></div>
<pre><code>## $ModelFit
##            lnL n.params       AIC       AICc
## shift.1 3.6472        3 -1.294399 -0.2943992
## 
## $Rates
##   node shiftPos          MLRate          LowerCI          UpperCI
## 1   45    clade 8.6033966785759 1.71610149858682 156.049851856678
## 
## $optimalTree
## 
## Phylogenetic tree with 28 tips and 27 internal nodes.
## 
## Tip labels:
##  A_alutaceu, A_inexpect, A_vanidicu, A_alfaroi, A_macilent, A_clivicol, ...
## Node labels:
##  2, 2, 2, 2, 2, 2, ...
## 
## Rooted; includes branch lengths.</code></pre>
<div class="sourceCode"><pre class="sourceCode r"><code class="sourceCode r">colour_motmot &lt;-<span class="st"> </span><span class="kw">plotPhylo.motmot</span>(<span class="dt">phy=</span>phy.clade, <span class="dt">traitMedusaObject=</span>trait.medusa.tm2.summary, <span class="dt">reconType =</span> <span class="st">&quot;rates&quot;</span>, <span class="dt">type =</span> <span class="st">&quot;fan&quot;</span>, <span class="dt">cex=</span><span class="fl">0.5</span>, <span class="dt">edge.width=</span><span class="dv">2</span>)</code></pre></div>
<div class="figure">
<img src="https://github.com/ghthomas/motmot/blob/motmot.2.0/vignettes/figures/plot9-1.png" alt="the subset of the tree" />
<p class="caption">the subset of the tree</p>
</div>
<p>Thomas and Freckleton (2012) showed the tm2 algortihm has a high type-one error rate. One way to ameriolate this is to estimate the level a one shift is supported when we know BM is the true model. For example, we could simulate 1000 BM datasets on the tree, estimate a single shift using the tm2 algortihm, and calculating the difference between the AICcs for each BM and one shift model. We can these use this difference to estimate the AICc ‘penalty’ the is needed to reduce the tm2 type-one error rate to 0.05. We could use this penalty in the ‘cutoff’ argument of the traitMedusaSummary argument.</p>
<p>This is shown but not run in the code below</p>
<div class="sourceCode"><pre class="sourceCode r"><code class="sourceCode r"><span class="co"># not run</span>
<span class="co"># sim.bm &lt;- transformPhylo.sim(phy=phy.clade, n=1000, model=&quot;bm&quot;)</span>
<span class="co"># aic.cut.off &lt;- apply(sim.bm, 2, function(x) {</span>
    <span class="co"># bm.test &lt;- transformPhylo.ML(y=as.matrix(x), phy=phy.clade, model=&quot;tm2&quot;, minCladeSize=2, nSplits=1)</span>
    <span class="co"># bm.test[[1]][,&quot;AICc&quot;]</span>
    <span class="co"># })</span>
<span class="co"># tm2.cut.off &lt;- quantile(aic.cut.off[1,] - aic.cut.off[2,], 0.95)</span></code></pre></div>
</div>
</div>
<div id="time-slice-model" class="section level1">
<h1>Time-slice model</h1>
<p>A new addition to motmot is a Maximum likelihood model that allows for heterogeneous rates in different times of evolution. These models are seperate from the models that allow for heterogeneous rates among lineages, as modelled by the ‘traitMedusa’ algorithms.</p>
<p>The ‘timeSlice’ model is implemented using the ‘transformPhylo.ML’ function, using the argument model = ‘timeSlice’. The function allows for two seperate models of evolution. In one, it is possible to test shifts in evolution at times selected <em>a priori</em>. Alternatively, the fit of models can be tested at a range of different times, and the function will return the best-fitting model</p>
<p>First we will test for a shift in the rate of evolution 10 million years ago.</p>
<div class="sourceCode"><pre class="sourceCode r"><code class="sourceCode r">timeSlice.<span class="fl">10.</span>ml &lt;-<span class="st"> </span><span class="kw">transformPhylo.ML</span>(<span class="dt">y=</span>male.length.clade, <span class="dt">phy=</span>phy.clade, <span class="dt">model=</span><span class="st">&quot;timeSlice&quot;</span>, <span class="dt">splitTime=</span><span class="kw">c</span>(<span class="dv">10</span>))</code></pre></div>
<pre><code>## [1] &quot;BM model&quot;
##         lnL         AIC        AICc    sigma.sq   anc.state 
## -0.28385540  4.56771081  5.04771081  0.00185807  3.84948140 
## [1] &quot;shiftModel&quot;
##          lnL          AIC         AICc     sigma.sq    anc.state 
##  2.946487446  2.107025107  3.846155542  0.001006388  3.860015287 
##        rate1        rate2  time.split1 
##  0.692072848  2.944767886 10.000000000</code></pre>
<p>We can use the function ‘timeSliceSummary’ to plot and summarise the results. The output summarises the best model according to AICc values. This function automatically plots the original tree showing the location of shift(s), and the colours show the relative rates in each time slice. The second plot below shows the same tree and colours, but with the branch lengths scaled to the ML optimised rates</p>
<div class="sourceCode"><pre class="sourceCode r"><code class="sourceCode r">outputSummary &lt;-<span class="st"> </span><span class="kw">timeSliceSummary</span>(timeSlice.<span class="fl">10.</span>ml, <span class="dt">cutoff=</span><span class="fl">0.001</span>, <span class="dt">cex.tip=</span><span class="fl">0.5</span>, <span class="dt">phylo.width=</span><span class="dv">2</span>, <span class="dt">colour.ramp=</span><span class="kw">c</span>(<span class="st">&quot;blue&quot;</span>, <span class="st">&quot;red&quot;</span>))</code></pre></div>
<div class="figure">
<img src="https://github.com/ghthomas/motmot/blob/motmot.2.0/vignettes/figures/plot10-1.png" alt="" />
<p class="caption">timeSlice plot</p>
</div>
<p>We can also see other summarise information, such as the CI for each rate estimate.</p>
<div class="sourceCode"><pre class="sourceCode r"><code class="sourceCode r">outputSummary<span class="op">$</span>RatesCI</code></pre></div>
<pre><code>##           rates       LCI       UCI
## rate1 0.6920728 0.1873338  2.115629
## rate2 2.9447679 0.9633096 10.877404</code></pre>
<p>Rather than testing the overall fit of each model, we can fit models to all times. The function automatically tests for all 1 Ma shifts between the age of the tree - 10 Ma, and the present + 10 Ma. We can specify a number of shifts we would like to test for. Here we will test for up to 3 shifts. The model will test one shift, save it, search for a second, save those two, etc…</p>
<p>Here will modify the boundary age argument so all split times are tested between 65-5 Myrs, using the ‘boundaryAge’ argument. As we are not tested set times we need to set the number of splits to test using ‘nSplits’ - we will allow up to 3 splits</p>
<div class="sourceCode"><pre class="sourceCode r"><code class="sourceCode r">timeSlice.ml &lt;-<span class="st"> </span><span class="kw">transformPhylo.ML</span>(<span class="dt">y=</span>male.length.clade, <span class="dt">phy=</span>phy.clade, <span class="dt">model=</span><span class="st">&quot;timeSlice&quot;</span>, <span class="dt">nSplits=</span><span class="dv">3</span>, <span class="dt">boundaryAge=</span><span class="dv">5</span>)</code></pre></div>
<pre><code>## [1] BM model
##         lnL         AIC        AICc    sigma.sq   anc.state 
## -0.28385540  4.56771081  5.04771081  0.00185807  3.84948140 
## [1] shift 1
##          lnL          AIC         AICc     sigma.sq    anc.state 
##  4.509949148 -1.019898297  0.719232138  0.001070818  3.858225598 
##       rates1       rates2  time.split1 
##  0.613247203  3.633771307  5.545642000 
## [1] shift 2
##           lnL           AIC          AICc      sigma.sq     anc.state 
##  7.7641301455 -5.5282602909 -2.8009875637  0.0006475698  3.8306699807 
##        rates1        rates2        rates3   time.split1   time.split2 
##  3.4703058871  0.0000000100  5.8887066401  5.5456420000 39.5456420000 
## [1] shift 3
##           lnL           AIC          AICc      sigma.sq     anc.state 
##  7.9487076095 -3.8974152190  0.1025847810  0.0006734181  3.8417362745 
##        rates1        rates2        rates3        rates4   time.split1 
##  0.0000000100  3.4550077401  0.0000000100  5.6540409593  5.5456420000 
##   time.split2   time.split3 
## 39.5456420000 65.5456420000</code></pre>
<p>And summarise the results. We can selected the cutoff AICc improvement needed to justify selecting the next model. Here we use the arbitary cut-off value of 1. We could test this formally by estimating the correct AICc value needed to reduced type-error &gt; 5% by using BM simulated data (an example using the tm2 is shown above)</p>
<div class="sourceCode"><pre class="sourceCode r"><code class="sourceCode r">outputSummary &lt;-<span class="st"> </span><span class="kw">timeSliceSummary</span>(timeSlice.ml, <span class="dt">cutoff=</span><span class="dv">1</span>, <span class="dt">cex.tip=</span><span class="fl">0.5</span>, <span class="dt">phylo.width=</span><span class="dv">2</span>, <span class="dt">colour.ramp=</span><span class="kw">c</span>(<span class="st">&quot;blue&quot;</span>, <span class="st">&quot;red&quot;</span>))</code></pre></div>
<div class="figure">
<img src="https://github.com/ghthomas/motmot/blob/motmot.2.0/vignettes/figures/plot11-1.png" alt="" />
<p class="caption">timeSlice plot</p>
</div>
<p>The output looks odd as the middle time period has be collapsed to zero</p>
</div>
<div id="nested-models-of-evolution" class="section level1">
<h1>Nested models of evolution</h1>
<p>We can also tested models of nested evolution in which an ancestral model of BM evolution changes to a alternative model (EB, OU, kappa, delta, psi) within the phylogeny (Puttick, In Review)</p>
<p>Here we can show an example of BM -&gt; OU and BM -&gt; ACDC at node 44 of the phylogeny. However, neither of these is significantly better than BM</p>
<div class="sourceCode"><pre class="sourceCode r"><code class="sourceCode r">bm.model &lt;-<span class="st"> </span><span class="kw">transformPhylo.ML</span>(male.length.clade, <span class="dt">phy=</span>phy.clade, <span class="dt">model=</span><span class="st">&quot;bm&quot;</span>)
nested.acdc &lt;-<span class="st"> </span><span class="kw">transformPhylo.ML</span>(male.length.clade, <span class="dt">phy=</span>phy.clade, <span class="dt">model=</span><span class="st">&quot;ACDC&quot;</span>, <span class="dt">nodeIDs=</span><span class="kw">c</span>(<span class="dv">44</span>))</code></pre></div>
<pre><code>## Warning in transformPhylo.ML(male.length.clade, phy = phy.clade, model =
## &quot;ACDC&quot;, : Confidence limits fall outside the current parameter bounds -
## consider changing lowerBound and/or upperBound</code></pre>
<div class="sourceCode"><pre class="sourceCode r"><code class="sourceCode r">nested.ou &lt;-<span class="st"> </span><span class="kw">transformPhylo.ML</span>(male.length.clade, <span class="dt">phy=</span>phy.clade, <span class="dt">model=</span><span class="st">&quot;OU&quot;</span>, <span class="dt">nodeIDs=</span><span class="kw">c</span>(<span class="dv">44</span>))</code></pre></div>
<pre><code>## Warning in transformPhylo.ML(male.length.clade, phy = phy.clade, model =
## &quot;OU&quot;, : Confidence limits fall outside parameter bounds - consider changing
## lowerBound and/or upperBound</code></pre>
<div class="sourceCode"><pre class="sourceCode r"><code class="sourceCode r"><span class="dv">1</span> <span class="op">-</span><span class="st"> </span><span class="kw">pchisq</span>(nested.acdc<span class="op">$</span>MaximumLikelihood <span class="op">-</span><span class="st"> </span>bm.model<span class="op">$</span>logLikelihood, <span class="dv">1</span>)</code></pre></div>
<pre><code>## [1] 0.05740847</code></pre>
<div class="sourceCode"><pre class="sourceCode r"><code class="sourceCode r"><span class="dv">1</span> <span class="op">-</span><span class="st"> </span><span class="kw">pchisq</span>(nested.ou<span class="op">$</span>MaximumLikelihood <span class="op">-</span><span class="st"> </span>bm.model<span class="op">$</span>logLikelihood, <span class="dv">1</span>)</code></pre></div>
<pre><code>## [1] 0.361424</code></pre>
</div>
<div id="bayesian-estimation-of-tree-transformation-models" class="section level1">
<h1>Bayesian estimation of tree transformation models</h1>
<p>The function ‘transformPhylo.MCMC’ allows for the estimation of model parameters using Bayesian statistics. Models of lambda, delta, kappa, OU, ACDC, and psi can currently be modelled using transformPhylo.MCMC</p>
<p>The model allows for a pre-optimisation step. The model we test 30 (default) different deviations for the acceptance proposal distribution in order for the model to achieve an acceptance of around 0.44. This is done by default in the model but can be turned off by setting ‘opt.accept.rate=FALSE’</p>
<p>We will run an MCMC chain of 2000 generations to estimate Pagel’s lambda and discarding the first 10% (‘200 generations (’burn.in = 0.1’). All the models use a ‘uniform’ prior for each of the parameters. For lambda, this is a uniform distribution between 0 and 1, meaning we think all potential values are equally likely</p>
<div class="sourceCode"><pre class="sourceCode r"><code class="sourceCode r"><span class="kw">set.seed</span>(<span class="dv">20</span>) <span class="co"># set seed so run will be identical - for example use only</span>
lambda.mcmc &lt;-<span class="st"> </span><span class="kw">transformPhylo.MCMC</span>(<span class="dt">y=</span>male.length.clade, <span class="dt">phy=</span>phy.clade, <span class="dt">model=</span><span class="st">&quot;lambda&quot;</span>, <span class="dt">mcmc.iteration=</span><span class="dv">2000</span>, <span class="dt">burn.in=</span><span class="fl">0.1</span>)</code></pre></div>
<pre><code>## optimising acceptance ratio fine-tune
##   running
 acceptance attempt 0.463 best acceptance 0 best SD 0.448
##  finished fine.tune
##   
 MCMC progress: 100.0000 %
## $median
##    Lambda 
## 0.7731313 
## 
## $`95.HPD`
## lower 95% HPD upper 95% HPD 
##     0.4845588     0.9449442 
## 
## $ESS
##   Lambda 
## 494.6253 
## 
## $acceptance.rate
## [1] 0.4708495</code></pre>
<p>We can know check the posterior estimate of lambda and convergence of the model. The median and 95 Highest Posterior Density (HPD) is output by the model. Some convergene diagnostics are output as standard: Effective Sample Size (ESS) and acceptance rate. We aim for an ESS of at least 200 and an acceptance rate around 0.44</p>
<div class="sourceCode"><pre class="sourceCode r"><code class="sourceCode r">lambda.mcmc[<span class="dv">1</span><span class="op">:</span><span class="dv">4</span>]</code></pre></div>
<pre><code>## $median
##    Lambda 
## 0.7731313 
## 
## $`95.HPD`
## lower 95% HPD upper 95% HPD 
##     0.4845588     0.9449442 
## 
## $ESS
##   Lambda 
## 494.6253 
## 
## $acceptance.rate
## [1] 0.4708495</code></pre>
<p>Our lambda median value is 0.77 but there is a large 95% HPD (0.48-0.95). The ESS and acceptance rate look ok. We can also plot the trace from the MCMC chain</p>
<div class="sourceCode"><pre class="sourceCode r"><code class="sourceCode r"><span class="kw">plot</span>(lambda.mcmc<span class="op">$</span>mcmc.chain, <span class="dt">type=</span><span class="st">&quot;l&quot;</span>, <span class="dt">ylim=</span><span class="kw">c</span>(<span class="dv">0</span>, <span class="dv">1</span>), <span class="dt">xlab=</span><span class="st">&quot;generations&quot;</span>, <span class="dt">ylab=</span><span class="st">&quot;lambda&quot;</span>, <span class="dt">las=</span><span class="dv">1</span>)</code></pre></div>
<div class="figure">
<img src="https://github.com/ghthomas/motmot/blob/motmot.2.0/vignettes/figures/plot12-1.png" alt="" />
<p class="caption">MCMC trace for Pagel’s lambda</p>
</div>
<p>We would hope this is flatish trace around the stable posterior distribution - it looks ok (if not great)</p>
<p>References</p>
<div id="refs" class="references">
<div id="ref-Blomberg2003">
<p>Blomberg, Simon P, Theodore Garland Jr, and Anthony R Ives. 2003. “Testing for phylogenetic signal in comparative data: behavorial traits more labile.” <em>Evolution</em> 57 (4): 717–45. doi:<a href="https://doi.org/10.1111/j.0014-3820.2003.tb00285.x">10.1111/j.0014-3820.2003.tb00285.x</a>.</p>
</div>
<div id="ref-Harmon2010">
<p>Harmon, Luke J, Jonathan B Losos, T Jonathan Davies, Rosemary G Gillespie, John L Gittleman, W Bryan Jennings, Kenneth H Kozak, et al. 2010. “Early burts of body size and shape evolution are rare in comparative data.” <em>Evolution</em> 64 (8): 2385–96. doi:<a href="https://doi.org/10.1111/j.1558-5646.2010.01025.x">10.1111/j.1558-5646.2010.01025.x</a>.</p>
</div>
<div id="ref-OMeara2006">
<p>O’Meara, Brian C, Cécile Ané, Michael J Sanderson, and Peter C Wainwright. 2006. “Testing for different rates of continuous trait evolution using likelihood.” <em>Evolution</em> 60 (5): 922–33. doi:<a href="https://doi.org/10.1111/j.0014-3820.2006.tb01171.x">10.1111/j.0014-3820.2006.tb01171.x</a>.</p>
</div>
<div id="ref-Thomas2012">
<p>Thomas, Gavin H, and Robert P Freckleton. 2012. “MOTMOT: Models of trait macroevolution on trees.” <em>Methods in Ecology and Evolution</em> 3 (1): 145–51. doi:<a href="https://doi.org/10.1111/j.2041-210X.2011.00132.x">10.1111/j.2041-210X.2011.00132.x</a>.</p>
</div>
<div id="ref-Thomas2006">
<p>Thomas, Gavin H, Robert P Freckleton, and Tamás Székely. 2006. “Comparative analyses of the influence of developmental mode on phenotypic diversification rates in shorebirds.” <em>Proceedings of the Royal Society B: Biological Sciences</em> 273 (1594): 1619–24. doi:<a href="https://doi.org/10.1098/rspb.2006.3488">10.1098/rspb.2006.3488</a>.</p>
</div>
</div>
</div>
</body>
</html>
