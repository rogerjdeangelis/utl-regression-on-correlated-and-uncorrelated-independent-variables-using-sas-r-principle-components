# utl-regression-on-correlated-and-uncorrelated-independent-variables-using-sas-r-principle-components
Regression on correlated and uncorrelated independent variables using sas r principle components
    %let pgm=utl-regression-on-correlated-and-uncorrelated-independent-variables-using-sas-r-principle-components;

    Regression on correlated and uncorrelated independent variables using sas r principle components

    github
    https://tinyurl.com/ynev6s4u
    https://github.com/rogerjdeangelis/utl-regression-on-correlated-and-uncorrelated-independent-variables-using-sas-r-principle-components

    SOAPBOX ON
    SAS proc reg is a work of mathematical art.
    Not surprising when the CEO of SAS is a Felllow of American Statistical Assoc
    and author of the sweep operator.
    SOAPBOX OFF

         Comparison with correlated independent variables vs uncorrelated regression

             1 Correlated independent variables X1 X2
             2 Uncorrelated independent variables PC1 PC2 (priciple components)
               Using principle components
               transform simutaneously X1 with X2 by applying PCA eigenvectors to X1 and X2
               Could also use Cholesky decomposition
             3 Perplexity summarization of Uncorrelated Regression on end
             4. Related repos

         The major drawback pf PCA regression is loss of information contained in the natural variables, ie age, height ...
         PCA creates two new variables wheree PC1 is a linear combination of X1 and X2 anf PC2 is an
         uncorrelated linear combination of X1 and X2.
         PCA accounts for the maximum variance in PC1 and lesser variance in PC2.

         The linear combination may have a natural interpretation. With a large number
         of variables PCA can reduce the dimensionality.

         Note the PCA transformations do not effect the predicted Y or the residuals Y-Yhat, but
         often reduces the variance of the betas.

         PCA regression tends to be more stable?

    /****************************************************************************************************************************/
    /*                                                            |                                                             */
    /*      CORRELATED X1 X2 DESIGN MATRIX                        |            UNCORRELATED PC1 PC2 DESIGN MATRIX               */
    /*      ==============================                        |            ==================================               */
    /*                                                            |                                                             */
    /*  SD1.HAVE total obs=30 23JUN2024:16:12:15                  |         SD1.WANT total obs=30 23JUN2024:16:13:25            */
    /*                                                            |                                                             */
    /*  Obs     Y       X1       X2                               |             Y          PC1         PC2                      */
    /*                                                            |                                                             */
    /*    1  1.24067  0.56680  0.15745                            |          1.24067    -1.16601    -0.32654                    */
    /*    2  0.68965  0.31237 -0.02737                            |          0.68965    -0.09071    -0.62645                    */
    /*    3 -0.67680 -0.55713 -0.05263                            |         -0.67680     1.32812     0.60443                    */
    /*    4 -0.30090 -0.07526 -0.15808                            |         -0.30090     0.98621    -0.52213                    */
    /*    5 -1.14115 -0.49168 -0.17685                            |         -1.14115     1.69057     0.04254                    */
    /*    6  0.29581  0.22145  0.11043                            |          0.29581    -0.46487     0.02475                    */
    /*                                                            |                                                             */
    /*--------------------------------------------------------------------------------------------------------------------------*/
    /*                                                            |                                                             */
    /*          CORRELATED X1 X2                                  |               UNCORRELATED PC1 PC2                          */
    /*          ----------------                                  |                                                             */
    /*  VARIABLE       X1         X2                              |          VARIABLE  PC1         PC2                          */
    /*                                                            |                                                             */
    /*     X1       1.00000     0.87812                           |           PC1   1.00000     0.00000                         */
    /*     X2       0.87812     1.00000                           |           PC2   0.00000     1.00000                         */
    /*                                                            |                                                             */
    /*--------------------------------------------------------------------------------------------------------------------------*/
    /*                                                            |                                                             */
    /* UNLIKE THE UNCORRELATED SOLUTION, MANUAL                   |  EASY MANUAL CALCULATION OF REGRESSION COEFICIENTS(BETAS)   */
    /* CALCULATION OF THE COEFICIENTS IS VERY DIFFICULT.          |                                                             */
    /*                                                            |  Lets use what SAS gives us to easily hand                  */
    /* Coeficients =Inverse(XTranspose * X)*(XTranspose*Y)        |  calculate the Regression coeficients                       */
    /*                                                            |                                                             */
    /* PROC REG MODEL CROSSPRODUCTS X'X X'Y Y'Y (REARRANGED)      |  Coeficients =Inverse(XTranspose * X)*(XTranspose*Y)        */
    /*             Xtranspose*X                XtransposeY        |                                                             */
    /*      ----------------------------------  ---------         |  From proc Reg Model Crossproducts Inv(X'X)*X'Y             */
    /* Var          X1           X2  Intercept    Y               |                  Xtranspose*X             XtransposeY       */
    /* X1    6.5577466 2.4132438137  3.0609460  12.0093350        |             ----------------------------  -------------     */
    /* X2    2.4132438 1.1370509087  1.6389147   4.7967337        |  Variable   Intercept   PC1           PC2             Y     */
    /* Inte  3.0609460 1.6389147656         30   6.9573554        |  Incpt      30            0             0  6.9573554362     */
    /*                                                            |  PC1         0 54.465399436             0  -33.64935772     */
    /* USING R                                                    |  PC2         0            0  3.5346005644  -0.784889585     */
    /* Xtranppose * X SAME AS PROC REG                            |                                                             */
    /*          X1       X2                                       |              XTranspose*X (from above)                      */
    /* X1 6.557747 2.413244  3.060946                             |  Icp 30              0                 0                    */
    /* X2 2.413244 1.137051  1.638915                             |  PC1  0   54.465399436                 0                    */
    /*    3.060946 1.638915 30.000000                             |  PC2  0             0       3.5346005644                    */
    /*                                                            |                                                             */
    /* inverse(Xtranspose * X)                                    |  Inverse is just the reciprocal of the diagonal elements.   */
    /*                         X1          X2                     |  The product of matrix * inverse is an identity matrix      */
    /* X1  0.69947527 -1.49977388  0.01056485                     |                                                             */
    /* X2 -1.49977388  4.17036644 -0.07480494                     |  Icp   1/30                     0                    0      */
    /*     0.01056485 -0.07480494  0.03634202                     |  PC1      0        1/54.465399436                    0      */
    /*                                                            |  PC2      0                     0       1/3.5346005644      */
    /* Xtranspose * Y (SAME AS PROC REG)                          |                                                      0      */
    /* X1 12.009335                                               |  The (XTranspose*Y from Proc Reg)                           */
    /* X2  4.796734                                               |  Incp   6.9573554362                                        */
    /*     6.957355                                               |  PC1    -33.64935772                                        */
    /*                                                            |  PC2    -0.784889585                                        */
    /* Coeficients=Inverse(XTranspose X)*(XTranspose*Y)           |                                                             */
    /* X1 1.27972030                                              |  Inverse of (XtransposeX) and XTransposeY                   */
    /* X2 1.47240586                                              |  Icp  0.03333333             0             0  6.9573554362  */
    /*    0.02090177                                              |  PC1           0  0.2239441318             0  -33.64935772  */
    /* Yhat = 1.27972030 + 1.47240586*X1 + 0.02090177*X2          |  PC2           0             0  0.2829173995  -0.784889585  */
    /*                                                            |                                                             */
    /*                                                            |  Coeficients= Inverse(XTranspose * X) * (XTranspose * Y )   */
    /*                                                            |                                                             */
    /*                                                            |  Intercept =  0.0333333333*6.9573554362 = 0.2319118456      */
    /*                                                            |       PC1  =  0.0183602803*-33.64935772 = -0.61781164       */
    /*                                                            |       PC2  =  0.2829173995*-0.784889585 =  -0.22205892      */
    /*                                                            |                                                             */
    /*                                                            |   Yhat = 0.2319118456 + -0.61781164 * PC1 -0.22205892*PC2   */
    /*                                                            |                                                             */
    /*--------------------------------------------------------------------------------------------------------------------------*/
    /*                                                            |                                                             */
    /*  Partitioning of SumSq depends on order                    |  Partitioning of SumSq does not depend on order             */
    /*  variables enter model                                     |  variables enter model                                      */
    /*  ---------------------------------------                   |  Very effective for selection techniques stepwise           */
    /*                                                            |  ------------------------------------------------           */
    /*                                                            |                                                             */
    /*                                                            |                                                             */
    /*          Results           Results                         |                                                             */
    /*     Effected by order   Not Efected(Adj)                   |      Uneffected by order   No need to adjust?               */
    /*     Type I SS           Type II SS                         |       Type I SS              Type II SS                     */
    /*                                                            |                                                             */
    /*  X1   20.44340              2.34130                        |    PC1   20.78896             20.78896                      */
    /*  X2   0.51985               0.51985                        |    PC2    0.17429             0.17429                       */
    /*                                                            |                                                             */
    /*------------------------------------------------------------------------------------------------------------------------  */
    /*                                                            |                                                             */
    /*  Tolerance                                                 |    Tolerance  Let p =%  of variace shared PC1 PC2)          */
    /*  ---------                                                 |    ---------      If p=0 then no shared variace and         */
    /*                                                            |                   tolerance=1                               */
    /*                                                            |                                                             */
    /*    0.22891  X1 modrate linear comb of X2                   |       1.00000       1 implies %variace shared PC1 PC2=0     */
    /*                                                            |    no colinearity      then muti-colinearity                */
    /*                                                            |                     0 -> PC1 is linear combination of PC2   */
    /*                                                            |                     0 -> PC1 is linear combination of PC2   */
    /*                                                            |                     0 -> PC1 is linear combination of PC2   */
    /*                                                            |                                                             */
    /*-----------------------------------------------------------------------------------------------------------------------   */
    /*                                                            |                                                             */
    /*   Variance                                                 |     Variance                                                */
    /*   Inflation                                                |    Inflation                                                */
    /*   ---------                                                |    ---------                                                */
    /*    4.36853 (1/tolerance)=1/0.22891                         |     1.00000         No multicollinearity                    */
    /*   VIF = 1: No multicollinearity                            |                                                             */
    /*   1 < VIF < 5: Moderate multicollinearity                  |                                                             */
    /*   VIF > 5: High multicollinearity                          |                                                             */
    /*    (may lead to unreliable beta estimates)                 |                                                             */
    /*                                                            |                                                             */
    /*------------------------------------------------------------------------------------------------------------------------  */
    /*                                                            |                                                             */
    /*  Wider standard errors                                     |      Tighter standars error (tighter confidence errors)     */
    /*                                                            |      --------------------------------------------------     */
    /*                Parameter  Standard                         |                    Parameter  Standard                      */
    /*  Variable  DF   Estimate    Error t Value                  |      Variable  DF   Estimate     Error  t Value             */
    /*                                                            |                                                             */
    /*  Intercept  1    0.02090  0.04453    0.47                  |      Intercept  1    0.23191   0.04264     5.44             */
    /*  X1         1    1.27972  0.19535    6.55                  |      PC1        1   -0.61781   0.03165   -19.52             */
    /*  X2         1    1.47241  0.47699    3.09                  |      PC2        1   -0.22206   0.12424    -1.79             */
    /*                                                            |                                                             */
    /* ---------------------------------------------------        |-------------------------------------------------------------*/
    /*                                                            |                                                             */
    /* SAS PROC REG OUTPUT                                        |        SAS PROC REG OUTPUT                                  */
    /* ===================                                        |        ===================                                  */
    /*                                                            |                                                             */
    /* Model Crossproducts X'X X'Y Y'Y                            |  Model Crossproducts X'X X'Y Y'Y   */                       */
    /*                                                            |                   Xtranspose*X               XtransposeY    */
    /* VAR    INTERCEPT            X1            X2            Y  |         ----------------------------------  ------------    */
    /*                Xtranspose*X                   XtransposeY  |  Var    Intercept          PC1          PC2            Y    */
    /*      --------------------------------------- ------------  |  Intercept    30             0           0  6.9573554362    */
    /* Int           30  3.0609460232  1.6389147656 6.9573554362  |  PC1            0 54.465399436              -33.64935772    */
    /* X1  3.0609460232  6.5577466614  2.4132438137 12.009335048  |  PC2            0            0 3.5346005644 -0.784889585    */
    /* X2  1.6389147656  2.4132438137  1.1370509087 4.7967337351  |  Y   6.9573554362 -33.64935772 -0.784889585 24.049797711    */
    /* Y   6.9573554362  12.009335048  4.7967337351 24.049797711  |                                                  |          */
    /*                                                            |                                                             */
    /* Dependent Variable: Y                                      |  Dependent Variable: Y                                      */
    /* Number of Observations Read          30                    |  Number of Observations Read          30                    */
    /*                                                            |                                                             */
    /*                 Analysis of Variance                       |                  Analysis of Variance                       */
    /*                         Sum of      Mean                   |                          Sum of     Mean                    */
    /* Source            DF   Squares    SquareF Value            |  Source            DF   Squares   Square F Value            */
    /*                                                            |                                                             */
    /* Model              2  20.96326  10.48163 192.12            |  Model              2  20.96326 10.48163  192.12            */
    /* Error             27   1.47305   0.05456                   |  Error             27   1.47305  0.05456                    */
    /* Corrected Total   29  22.43630                             |  Corrected Total   29  22.43630                             */
    /*                                                            |                                                             */
    /*    Analysis of Variance                                    |     Analysis of Variance                                    */
    /* Source              Pr > F                                 |  Source              Pr > F                                 */
    /*                                                            |                                                             */
    /* Model               <.0001                                 |  Model               <.0001                                 */
    /* Error                                                      |  Error                                                      */
    /* Corrected Total                                            |  Corrected Total                                            */
    /*                                                            |                                                             */
    /*                                                            |  Root MSE              0.23358    R-Square  0.9343          */
    /* Root MSE              0.23358    R-Square  0.9343          |  Dependent Mean        0.23191    Adj R-Sq  0.9295          */
    /* Dependent Mean        0.23191    Adj R-Sq  0.9295          |  Coeff Var           100.71720                              */
    /* Coeff Var           100.71720                              |                                                             */
    /*                                                            |                                                             */
    /*                  Parameter Estimates                       |                        Parameter Estimates                  */
    /*               Parameter  Standard                          |                Parameter  Standard                          */
    /* Variable  DF   Estimate    Error t Value Pr > |t|          |  Variable  DF   Estimate     Error  t Value Pr > |t|        */
    /*                                                            |                                                             */
    /* Intercept  1    0.02090  0.04453    0.47   0.6425          |  Intercept  1    0.23191   0.04264     5.44   <.0001        */
    /* X1         1    1.27972  0.19535    6.55   <.0001          |  PC1        1   -0.61781   0.03165   -19.52   <.0001        */
    /* X2         1    1.47241  0.47699    3.09   0.0046          |  PC2        1   -0.22206   0.12424    -1.79   0.0851        */
    /*                                                            |                                                             */
    /*               Parameter Estimates                          |                     Parameter Estimates                     */
    /* Variable  DF  Type I SS    Type II SS  Tolerance           |  Variable    DF     Type I SS Type II SS Tolerance          */
    /*                                                            |                                                             */
    /* Intercept  1    1.61349    0.01202           .             |  Intercept    1       1.61349    1.61349         .          */
    /* X1         1   20.44340    2.34130     0.22891             |  PC1          1      20.78896   20.78896   1.00000          */
    /* X2         1    0.51985    0.51985     0.22891             |  PC2          1       0.17429    0.17429   1.00000          */
    /*                                                            |                                                             */
    /*     Parameter Estimates                                    |      Parameter Estimates                                    */
    /*                     Variance                               |                      Variance                               */
    /* Variable    DF     Inflation                               |  Variable    DF     Inflation                               */
    /*                                                            |                                                             */
    /* Intercept    1             0                               |  Intercept    1             0                               */
    /* X1           1       4.36853                               |  PC1          1       1.00000                               */
    /* X2           1       4.36853                               |  PC2          1       1.00000                               */
    /*                                                            |                                                             */
    /*                                                            |                                                             */
    /*         Correlation of Estimates                           |        Correlation of Estimates                             */
    /* Variable        X1                X2                       |                                                             */
    /*                                                            |  Variable       PC1        PC2                              */
    /* X1          1.0000           -0.8781                       |                                                             */
    /* X2         -0.8781            1.0000                       |  PC1         1.0000     0.0000                              */
    /*                                                            |                                                             */
    /****************************************************************************************************************************/

    /*                           _       _           _                                   _
    / |   ___ ___  _ __ _ __ ___| | __ _| |_ ___  __| | _ __ ___  __ _ _ __ ___  ___ ___(_) ___  _ __
    | |  / __/ _ \| `__| `__/ _ \ |/ _` | __/ _ \/ _` || `__/ _ \/ _` | `__/ _ \/ __/ __| |/ _ \| `_ \
    | | | (_| (_) | |  | | |  __/ | (_| | ||  __/ (_| || | |  __/ (_| | | |  __/\__ \__ \ | (_) | | | |
    |_|  \___\___/|_|  |_|  \___|_|\__,_|\__\___|\__,_||_|  \___|\__, |_|  \___||___/___/_|\___/|_| |_|
     _                   _                                       |___/
    (_)_ __  _ __  _   _| |_
    | | `_ \| `_ \| | | | __|
    | | | | | |_) | |_| | |_
    |_|_| |_| .__/ \__,_|\__|
            |_|
    */

    options validvarname=upcase;
    libname sd1 "d:/sd1";
    data sd1.have;
      call streaminit(4321);
      do i=1 to 30;
        y=0.05 + rand('normal',0,1);
        x1 = y * 0.5 + rand('normal',0,0.1);
        x2 = y * 0.3 - x1 * 0.2 + rand('normal',0,0.1);
        output;
      end;
      drop i;
    run;quit;

    /**************************************************************************************************************************/
    /*                                                                                                                        */
    /*  SD1.HAVE total obs=30                                                                                                 */
    /*                                                                                                                        */
    /*  Obs        Y          X1          X2                                                                                  */
    /*                                                                                                                        */
    /*    1     1.29067     0.59180     0.16745                                                                               */
    /*    2     0.73965     0.33737    -0.01737                                                                               */
    /*    3    -0.62680    -0.53213    -0.04263                                                                               */
    /*    ...                                                                                                                 */
    /*   28    -0.74911    -0.41607    -0.10870                                                                               */
    /*   29    -1.02584    -0.34066    -0.22165                                                                               */
    /*   30     1.82865     0.76747     0.46834                                                                               */
    /*                                                                                                                        */
    /**************************************************************************************************************************/

    /*
     _ __  _ __ ___   ___ ___  ___ ___
    | `_ \| `__/ _ \ / __/ _ \/ __/ __|
    | |_) | | | (_) | (_|  __/\__ \__ \
    | .__/|_|  \___/ \___\___||___/___/
    |_|
    */

    options ls=120 ps=255;
    proc reg data=sd1.have;
     Correlated: model y=x1 x2 /vif tol corrb covb  gmsep ss1 ss2 xpx;
     output out=Pred pred=yhat;
    run;quit;

    /*----                                                                   ----*/
    /*---- Use R to show how to get the coef using just proc reg output      ----*/
    /*----                                                                   ----*/

    %utl_submit_r64x('
    library(haven);
    source("c:/oto/fn_tosas9x.R");
    have<-as.matrix(read_sas("d:/sd1/have.sas7bdat"));
    have;
    X<-have[,c(2,3)];
    X <- cbind(X, rep(1, nrow(X)));
    Y<-have[,1];
    Y;
    XtX <- t(X) %*% X;
    print("X transpose * X (SAME AS SAS PROC REG");
    XtX;
    XtY <- t(X) %*% Y;
    print("X transpose * Y (SAME AS SAS PROC REG");
    XtY;
    print("inverse (X transpose * X)");
    XtX_inv <- solve(XtX);
    XtX_inv;
    print("Betas as Inverse(XTranspose X)*(XTranspose*Y");
    b <- XtX_inv %*% XtY;
    print(b);
    ');

    /**************************************************************************************************************************/
    /*  _ __  _ __ ___   ___   _ __ ___  __ _                                                                                 */
    /* | `_ \| `__/ _ \ / __| | `__/ _ \/ _` |                                                                                */
    /* | |_) | | | (_) | (__  | | |  __/ (_| |                                                                                */
    /* | .__/|_|  \___/ \___| |_|  \___|\__, |                                                                                */
    /* |_|                              |___/                                                                                 */
    /*                                                                                                                        */
    /* Model: Correlated                                                                                                      */
    /*                                                                                                                        */
    /*                          Model Crossproducts X'X X'Y Y'Y                                                               */
    /*                                                                                                                        */
    /* Variable          Intercept                X1                X2                 Y                                      */
    /*                                                                                                                        */
    /* Intercept                30      3.0609460232      1.6389147656      6.9573554362                                      */
    /* X1             3.0609460232      6.5577466614      2.4132438137      12.009335048                                      */
    /* X2             1.6389147656      2.4132438137      1.1370509087      4.7967337351                                      */
    /* Y              6.9573554362      12.009335048      4.7967337351      24.049797711                                      */
    /*                                                                                                                        */
    /*                                                                                                                        */
    /* The REG Procedure                                                                                                      */
    /* Model: Correlated                                                                                                      */
    /* Dependent Variable: Y                                                                                                  */
    /*                                                                                                                        */
    /* Number of Observations Read          30                                                                                */
    /* Number of Observations Used          30                                                                                */
    /*                                                                                                                        */
    /*                              Analysis of Variance                                                                      */
    /*                                     Sum of           Mean                                                              */
    /* Source                   DF        Squares         Square    F Value    Pr > F                                         */
    /*                                                                                                                        */
    /* Model                     2       20.96326       10.48163     192.12    <.0001                                         */
    /* Error                    27        1.47305        0.05456                                                              */
    /* Corrected Total          29       22.43630                                                                             */
    /*                                                                                                                        */
    /*                                                                                                                        */
    /* Root MSE              0.23358    R-Square     0.9343                                                                   */
    /* Dependent Mean        0.23191    Adj R-Sq     0.9295                                                                   */
    /* Coeff Var           100.71720                                                                                          */
    /*                                                                                                                        */
    /*                                                   Parameter Estimates                                                  */
    /*              Parameter Standard                                                  Variance                              */
    /* Variable  DF  Estimate    Error t Value Pr > |t| Type I SS Type II SS Tolerance  Inflation                             */
    /*                                                                                                                        */
    /* Intercept  1   0.02090  0.04453    0.47   0.6425   1.61349    0.01202         .          0                             */
    /* X1         1   1.27972  0.19535    6.55   <.0001  20.44340    2.34130   0.22891    4.36853                             */
    /* X2         1   1.47241  0.47699    3.09   0.0046   0.51985    0.51985   0.22891    4.36853                             */
    /*                                                                                                                        */
    /*                                                                                                                        */
    /*                              Covariance of Estimates                                                                   */
    /* Variable          Intercept                X1                X2                                                        */
    /*                                                                                                                        */
    /* Intercept      0.0019827234      0.0005763901      -0.004081158                                                        */
    /* X1             0.0005763901      0.0381615043      -0.081823661                                                        */
    /* X2             -0.004081158      -0.081823661      0.2275240664                                                        */
    /*                                                                                                                        */
    /*                                                                                                                        */
    /*                    Correlation of Estimates                                                                            */
    /* Variable          Intercept                X1                X2                                                        */
    /*                                                                                                                        */
    /* Intercept            1.0000            0.0663           -0.1921                                                        */
    /* X1                   0.0663            1.0000           -0.8781                                                        */
    /* X2                  -0.1921           -0.8781            1.0000                                                        */
    /*                                                                                                                        */
    /*------------------------------------------------------------------------------------------------------------------------*/
    /*  ____                                                                                                                  */
    /* |  _ \                                                                                                                 */
    /* | |_) |                                                                                                                */
    /* |  _ <                                                                                                                 */
    /* |_| \_\                                                                                                                */
    /*                                                                                                                        */
    /*  [1] "X transpose * X (SAME AS SAS PROC REG)"                                                                          */
    /*           X1       X2                                                                                                  */
    /*  X1 6.557747 2.413244  3.060946                                                                                        */
    /*  X2 2.413244 1.137051  1.638915                                                                                        */
    /*     3.060946 1.638915 30.000000                                                                                        */
    /*                                                                                                                        */
    /*  [1] "X transpose * Y (SAME AS SAS PROC REG)"                                                                          */
    /*          [,1]                                                                                                          */
    /*  X1 12.009335                                                                                                          */
    /*  X2  4.796734                                                                                                          */
    /*      6.957355                                                                                                          */
    /*                                                                                                                        */
    /*  [1] "inverse (X transpose * X)"                                                                                       */
    /*              X1          X2                                                                                            */
    /*  X1  0.69947527 -1.49977388  0.01056485                                                                                */
    /*  X2 -1.49977388  4.17036644 -0.07480494                                                                                */
    /*      0.01056485 -0.07480494  0.03634202                                                                                */
    /*                                                                                                                        */
    /*  [1] "Betas as Inverse(XTranspose X)*(XTranspose*Y"                                                                    */
    /*           [,1]                                                                                                         */
    /*  X1 1.27972030                                                                                                         */
    /*  X2 1.47240586                                                                                                         */
    /*     0.02090177                                                                                                         */
    /*                                                                                                                        */
    /*                                                                                                                        */
    /*  Yhat = 1.27972030 + 1.47240586*X1 + 0.02090177*X2                                                                     */
    /*                                                                                                                        */
    /**************************************************************************************************************************/

    /*___                                        _       _           _                                   _
    |___ \   _   _ _ __   ___ ___  _ __ _ __ ___| | __ _| |_ ___  __| | _ __ ___  __ _ _ __ ___  ___ ___(_) ___  _ __
      __) | | | | | `_ \ / __/ _ \| `__| `__/ _ \ |/ _` | __/ _ \/ _` || `__/ _ \/ _` | `__/ _ \/ __/ __| |/ _ \| `_ \
     / __/  | |_| | | | | (_| (_) | |  | | |  __/ | (_| | ||  __/ (_| || | |  __/ (_| | | |  __/\__ \__ \ | (_) | | | |
    |_____|  \__,_|_| |_|\___\___/|_|  |_|  \___|_|\__,_|\__\___|\__,_||_|  \___|\__, |_|  \___||___/___/_|\___/|_| |_|
     _                   _                                                         |___/
    (_)_ __  _ __  _   _| |_
    | | `_ \| `_ \| | | | __|
    | | | | | |_) | |_| | |_
    |_|_| |_| .__/ \__,_|\__|
            |_|
    */

    /*----                                                                   ----*/
    /*---- Transform original X1 and X2 to uncorrelated PC1 and PC2          ----*/
    /*---- PC1 and PC2 are uncorrelated linerar combinations of X1 and X2    ----*/
    /*---- PC1 accounts for the maximum variance in the X dsfign matrix      ----*/
    /*----                                                                   ----*/

    %utl_submit_r64x('
    library(stats);
    library(haven);
    source("c:/oto/fn_tosas9x.R");
    have<-read_sas("d:/sd1/have.sas7bdat");
    head(have);
    data <- data.frame(have$X1, have$X2);
    cor(data);
    pca <- prcomp(data, center = TRUE, scale. = TRUE);
    uncor <- as.data.frame(pca$x);
    cor(uncor);
    addy<-cbind(Y=have$Y,uncor);
    fn_tosas9x(
          inp    = addy
         ,outlib ="d:/sd1/"
         ,outdsn ="want"
         );
    ');

    proc print data=sd1.want;
    run;quit;

    /**************************************************************************************************************************/
    /*                                                                                                                        */
    /* Obs    ROWNAMES        Y          PC1         PC2                                                                      */
    /*                                                                                                                        */
    /*   1        1        1.29067    -1.16601    -0.32654                                                                    */
    /*   2        2        0.73965    -0.09071    -0.62645                                                                    */
    /*   3        3       -0.62680     1.32812     0.60443                                                                    */
    /*  ...                                                                                                                   */
    /*  28       28       -0.74911     1.39713     0.18175                                                                    */
    /*  29       29       -1.02584     1.70244    -0.35338                                                                    */
    /*  30       30        1.82865    -2.55316     0.52527                                                                    */
    /*                                                                                                                        */
    /**************************************************************************************************************************/
    /*
     _ __  _ __ ___   ___ ___  ___ ___
    | `_ \| `__/ _ \ / __/ _ \/ __/ __|
    | |_) | | | (_) | (_|  __/\__ \__ \
    | .__/|_|  \___/ \___\___||___/___/
    |_|
    */

    options ls=110 ps=255;
    proc reg data=sd1.want;
     Uncorrelated: model y=pc1 pc2 /vif tol corrb covb  gmsep ss1 ss2 xpx;
     output out=pcPred pred=pcyhat;
    run;quit;

    /*----                                                                   ----*/
    /*---- Check with R                                                      ----*/
    /*---- Oneliner in (solve(t(X) %*% X) %* % t(X) %*% Y)                   ----*/
    /*----                                                                   ----*/

    %utl_submit_r64x('
    library(haven);
    source("c:/oto/fn_tosas9x.R");
    have<-as.matrix(read_sas("d:/sd1/want.sas7bdat"));
    have;
    X<-have[,c(3,4)];
    X <- cbind(X, rep(1, nrow(X)));
    Y<-have[,2];
    XtX <- t(X) %*% X;
    print("X transpose * X (SAME AS SAS PROC REG");
    XtX;
    XtY <- t(X) %*% Y;
    print("X transpose * Y (SAME AS SAS PROC REG");
    XtY;
    print("inverse (X transpose * X)");
    XtX_inv <- solve(XtX);
    XtX_inv;
    print("Betas as Inverse(XTranspose X)*(XTranspose*Y");
    b <- XtX_inv %*% XtY;
    print(b);

    ');
    /**************************************************************************************************************************/
    /*                                                                                                                        */
    /*  _ __  _ __ ___   ___   _ __ ___  __ _                                                                                 */
    /* | `_ \| `__/ _ \ / __| | `__/ _ \/ _` |                                                                                */
    /* | |_) | | | (_) | (__  | | |  __/ (_| |                                                                                */
    /* | .__/|_|  \___/ \___| |_|  \___|\__, |                                                                                */
    /* |_|                              |___/                                                                                 */
    /*                                                                                                                        */
    /* The REG Procedure                                                                                                      */
    /* Model: Uncorrelated                                                                                                    */
    /*                                                                                                                        */
    /*                      Model Crossproducts X'X X'Y Y'Y                                                                   */
    /*                                  X transpose X                  X transpose Y                                          */
    /*               --------------------------------------------      --------------                                         */
    /* Variable      Intercept               PC1               PC2                 Y                                          */
    /* Intercept            30                 0                0       6.9573554362                                          */
    /* PC1                   0      54.465399436                0       -33.64935772                                          */
    /* PC2                   0                 0      3.5346005644      -0.784889585                                          */
    /*                                                                                                                        */
    /* Y          6.9573554362      -33.64935772      -0.784889585      24.049797711                                          */
    /*                                                                                                                        */
    /* Dependent Variable: Y                                                                                                  */
    /*                                                                                                                        */
    /* Number of Observations Read          30                                                                                */
    /* Number of Observations Used          30                                                                                */
    /*                                                                                                                        */
    /*                                                                                                                        */
    /*                              Analysis of Variance                                                                      */
    /*                                                                                                                        */
    /*                                     Sum of           Mean                                                              */
    /* Source                   DF        Squares         Square    F Value    Pr > F                                         */
    /*                                                                                                                        */
    /* Model                     2       20.96326       10.48163     192.12    <.0001                                         */
    /* Error                    27        1.47305        0.05456                                                              */
    /* Corrected Total          29       22.43630                                                                             */
    /*                                                                                                                        */
    /*                                                                                                                        */
    /* Root MSE              0.23358    R-Square     0.9343                                                                   */
    /* Dependent Mean        0.23191    Adj R-Sq     0.9295                                                                   */
    /* Coeff Var           100.71720                                                                                          */
    /*                                                                                                                        */
    /*                                                                                                                        */
    /*                                              Parameter Estimates                                                       */
    /*                                                                                                                        */
    /*               Parameter Standard                                                  Variance                             */
    /* Variable   DF  Estimate    Error t Value Pr > |t| Type I SS Type II SS Tolerance Inflation                             */
    /*                                                                                                                        */
    /* Intercept   1   0.23191  0.04264    5.44   <.0001   1.61349    1.61349         .         0                             */
    /* PC1         1  -0.61781  0.03165  -19.52   <.0001  20.78896   20.78896   1.00000   1.00000                             */
    /* PC2         1  -0.22206  0.12424   -1.79   0.0851   0.17429    0.17429   1.00000   1.00000                             */
    /*                                                                                                                        */
    /*                                                                                                                        */
    /*                  Covariance of Estimates                                                                               */
    /*                                                                                                                        */
    /* Variable          Intercept               PC1               PC2                                                        */
    /*                                                                                                                        */
    /* Intercept      0.0018185777      4.448388E-20      3.427307E-19                                                        */
    /* PC1            4.448388E-20      0.0010016879      -4.02728E-18                                                        */
    /* PC2            3.427307E-19      -4.02728E-18      0.0154352185                                                        */
    /*                                                                                                                        */
    /*                                                                                                                        */
    /*                    Correlation of Estimates                                                                            */
    /*                                                                                                                        */
    /* Variable          Intercept               PC1               PC2                                                        */
    /*                                                                                                                        */
    /* Intercept            1.0000            0.0000            0.0000                                                        */
    /* PC1                  0.0000            1.0000           -0.0000                                                        */
    /* PC2                  0.0000           -0.0000            1.0000                                                        */
    /*                                                                                                                        */
    /*------------------------------------------------------------------------------------------------------------------------*/
    /*  ____                                                                                                                  */
    /* |  _ \                                                                                                                 */
    /* | |_) |                                                                                                                */
    /* |  _ <                                                                                                                 */
    /* |_| \_\                                                                                                                */
    /*                                                                                                                        */
    /* [1] "X transpose * X (SAME AS SAS PROC REG"                                                                            */
    /*               PC1           PC2                                                                                        */
    /* PC1  5.446540e+01             0             0                                                                          */
    /* PC2             )  3.534601e+00             0                                                                          */
    /*                 0             0  3.000000e+01                                                                          */
    /*                                                                                                                        */
    /* [1] "X transpose * Y (SAME AS SAS PROC REG"                                                                            */
    /*                                                                                                                        */
    /* PC1 -33.6493577                                                                                                        */
    /* PC2  -0.7848896                                                                                                        */
    /*       6.9573554                                                                                                        */
    /*                                                                                                                        */
    /* [1] "inverse (X transpose * X)"                                                                                        */
    /*                                                                                                                        */
    /*               PC1           PC2                                                                                        */
    /* PC1  1.836028e-02             0            0                                                                           */
    /* PC2             0  2.829174e-01            0                                                                           */
    /*                 0             0 3.333333e-02                                                                           */
    /*                                                                                                                        */
    /* [1] "Betas as Inverse(XTranspose X)*(XTranspose*Y"                                                                     */
    /*                                                                                                                        */
    /* PC1 -0.6178116                                                                                                         */
    /* PC2 -0.2220589                                                                                                         */
    /*      0.2319118                                                                                                         */
    /*                                                                                                                        */
    /**************************************************************************************************************************/

    /*____                         _           _ _
    |___ /   _ __   ___ _ __ _ __ | | _____  _(_) |_ _   _
      |_ \  | `_ \ / _ \ `__| `_ \| |/ _ \ \/ / | __| | | |
     ___) | | |_) |  __/ |  | |_) | |  __/>  <| | |_| |_| |
    |____/  | .__/ \___|_|  | .__/|_|\___/_/\_\_|\__|\__, |
            |_|             |_|                      |___/
    */
    When using uncorrelated independent variables in linear regression, several important
    changes occur compared to using correlated variables:

    1. Improved coefficient interpretation: With uncorrelated independent variables, each
     coefficient represents the unique effect of that variable on the dependent variable,
     holding all other variables constant. This makes interpretation more straightforward and reliable[1][4].

    2. Increased precision of coefficient estimates: Uncorrelated variables lead to more
     precise estimates of regression coefficients. The standard errors of the coefficients
     tend to be smaller, resulting in narrower confidence intervals[2].

    3. Reduced multicollinearity: By definition, uncorrelated variables eliminate
     multicollinearity issues. This avoids problems such as unstable coefficient
     estimates and inflated standard errors that can occur with correlated predictors[4].

    4. Simplified variable selection: When variables are uncorrelated, the inclusion
     or exclusion of one variable does not significantly affect the coefficients of
     other variables. This makes variable selection and model building more straightforward[1].

    5. Improved model stability: Uncorrelated variables lead to more stable models,
     as the coefficient estimates are less likely to change dramatically with small
     changes in the data or model specification[2].

    6. Orthogonal design: With uncorrelated variables, the design matrix becomes orthogonal.
     This simplifies many mathematical properties of the regression model and
     can lead to computational advantages[5].

    7. Easier assessment of variable importance: When variables are uncorrelated,
     their relative importance in predicting the dependent variable can be more
     directly compared using standardized coefficients (beta weights) or the increment
     in R-squared when each variable is added to the model[5].

    8. Reduced risk of overfitting: Uncorrelated variables are less likely to introduce
     redundant information into the model, potentially reducing the risk of overfitting,
     especially in smaller datasets[1].

    It's important to note that while using uncorrelated variables can simplify
     interpretation and improve certain aspects of the regression model, it may not
     always be possible or desirable in real-world scenarios. Many naturally occurring
     variables are correlated, and forcing orthogonality through data transformation or
     variable selection may lead to loss of important information or
    reduced predictive power.

    Citations:
    [1] https://stats.stackexchange.com/questions/455941/in-the-case-of-linear-regression-if-the-parameters-are-uncorrelated-does-this
    [2] https://timeseriesreasoning.com/contents/effect-of-irrelevant-variables/
    [3] https://statproofbook.github.io/P/slr-olscorr.html
    [4] https://statisticsbyjim.com/regression/multicollinearity-in-regression-analysis/
    [5] http://faculty.cas.usf.edu/mbrannick/regression/Part3/ImportanceNarrative.html
    /*  _                              
    | || |    _ __ ___ _ __   ___  ___ 
    | || |_  | `__/ _ \ `_ \ / _ \/ __|
    |__   _| | | |  __/ |_) | (_) \__ \
       |_|   |_|  \___| .__/ \___/|___/
                  |_|              
    */
    REPO                                                                                                                              
    ----------------------------------------------------------------------------------------------------------------------------------
    https://github.com/rogerjdeangelis/utl-betas-for-rolling-regressions                                                              
    https://github.com/rogerjdeangelis/utl-calculate-regression-coeficients-in-base-sas-fcmp-proc-reg-r-and-python                    
    https://github.com/rogerjdeangelis/utl-calculate-the-regression-slope-for-each-patient-by-treatment                               
    https://github.com/rogerjdeangelis/utl-drop-down-to-python-for-a-regression-sas-python-interface                                  
    https://github.com/rogerjdeangelis/utl-generate-all-possible-paiwise-interactions-products-regression                             
    https://github.com/rogerjdeangelis/utl-linear-regression-in-python-R-and-sas                                                      
    https://github.com/rogerjdeangelis/utl-locating-breakpoints-for-dogleg-mutiple-regression-lines                                   
    https://github.com/rogerjdeangelis/utl-outlier-analysis-based-on-robust-regression                                                
    https://github.com/rogerjdeangelis/utl-piecewise-regression-find-the-breakpoint                                                   
    https://github.com/rogerjdeangelis/utl-random-forest-regression-vs-linear-regression-with-uncorrelated-independent-variables-in-r 
    https://github.com/rogerjdeangelis/utl-regression-line-plus-and-minus-the-interquartile-range-of-dependent-variable               
    https://github.com/rogerjdeangelis/utl-scatter-plot-with-regression-line-coefficients-and-pvalue-in-one-datastep-sgplot           
    https://github.com/rogerjdeangelis/utl-simple-example-of-meta-regression-using-SAS-and-R                                          
    https://github.com/rogerjdeangelis/utl-using-linear-regression-with-base-sas-and-r-to-interpolate-missimg-values                  
    https://github.com/rogerjdeangelis/utl-using-the-regression-equation-to-score-a-table-like_a_weighted-sum                         
    https://github.com/rogerjdeangelis/utl_dosubl_do_regressions_when_data_is_between_dates                                           
    https://github.com/rogerjdeangelis/utl_excluding_rolling_regressions_with_one_on_more_missing_values_in_the_window                
    https://github.com/rogerjdeangelis/utl_how_to_automate_a_series_of_logistic_regressions                                           
    https://github.com/rogerjdeangelis/utl_multiple-regressions-using-arrays                                                          

    /*              _
      ___ _ __   __| |
     / _ \ `_ \ / _` |
    |  __/ | | | (_| |
     \___|_| |_|\__,_|

    */

