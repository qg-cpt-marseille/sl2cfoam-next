#include <math.h>
#include <quadmath.h>

#include "integration_gk.h"
#include "utils.h"
#include "error.h"

// Gauss quadrature weights and kronrod quadrature abscissae and
// weights as evaluated with 80 decimal digit arithmetic.
// Adapted from gsl source 
// (or see https://www.advanpix.com/2011/11/07/gauss-kronrod-quadrature-nodes-weights/)

// all values given in (0 1), must double to cover all (-1 1)
// xgk[1], xgk[3], ... abscissae of the N-point gauss rule. 
// xgk[0], xgk[2], ... abscissae of Kronrod extension
static const __float128 xgk_61[31] =
{
    0.999484410050490637571325895705811Q,
    0.996893484074649540271630050918695Q,
    0.991630996870404594858628366109486Q,
    0.983668123279747209970032581605663Q,
    0.973116322501126268374693868423707Q,
    0.960021864968307512216871025581798Q,
    0.944374444748559979415831324037439Q,
    0.926200047429274325879324277080474Q,
    0.905573307699907798546522558925958Q,
    0.882560535792052681543116462530226Q,
    0.857205233546061098958658510658944Q,
    0.829565762382768397442898119732502Q,
    0.799727835821839083013668942322683Q,
    0.767777432104826194917977340974503Q,
    0.733790062453226804726171131369528Q,
    0.697850494793315796932292388026640Q,
    0.660061064126626961370053668149271Q,
    0.620526182989242861140477556431189Q,
    0.579345235826361691756024932172540Q,
    0.536624148142019899264169793311073Q,
    0.492480467861778574993693061207709Q,
    0.447033769538089176780609900322854Q,
    0.400401254830394392535476211542661Q,
    0.352704725530878113471037207089374Q,
    0.304073202273625077372677107199257Q,
    0.254636926167889846439805129817805Q,
    0.204525116682309891438957671002025Q,
    0.153869913608583546963794672743256Q,
    0.102806937966737030147096751318001Q,
    0.051471842555317695833025213166723Q,
    0.000000000000000000000000000000000Q
};

// weights for Gauss points
static const double wg_61[15] =
{
    0.007968192496166605615465883474674,
    0.018466468311090959142302131912047,
    0.028784707883323369349719179611292,
    0.038799192569627049596801936446348,
    0.048402672830594052902938140422808,
    0.057493156217619066481721689402056,
    0.065974229882180495128128515115962,
    0.073755974737705206268243850022191,
    0.080755895229420215354694938460530,
    0.086899787201082979802387530715126,
    0.092122522237786128717632707087619,
    0.096368737174644259639468626351810,
    0.099593420586795267062780282103569,
    0.101762389748405504596428952168554,
    0.102852652893558840341285636705415
};

// weights for Kronrod points
static const double wgk_61[31] =
{
    0.001389013698677007624551591226760,
    0.003890461127099884051267201844516,
    0.006630703915931292173319826369750,
    0.009273279659517763428441146892024,
    0.011823015253496341742232898853251,
    0.014369729507045804812451432443580,
    0.016920889189053272627572289420322,
    0.019414141193942381173408951050128,
    0.021828035821609192297167485738339,
    0.024191162078080601365686370725232,
    0.026509954882333101610601709335075,
    0.028754048765041292843978785354334,
    0.030907257562387762472884252943092,
    0.032981447057483726031814191016854,
    0.034979338028060024137499670731468,
    0.036882364651821229223911065617136,
    0.038678945624727592950348651532281,
    0.040374538951535959111995279752468,
    0.041969810215164246147147541285970,
    0.043452539701356069316831728117073,
    0.044814800133162663192355551616723,
    0.046059238271006988116271735559374,
    0.047185546569299153945261478181099,
    0.048185861757087129140779492298305,
    0.049055434555029778887528165367238,
    0.049795683427074206357811569379942,
    0.050405921402782346840893085653585,
    0.050881795898749606492297473049805,
    0.051221547849258772170656282604944,
    0.051426128537459025933862879215781,
    0.051494729429451567558340433647099
};

// Compute the abscissae for the Gauss-Kronrod method.
// NB: the same pattern is then used in gk_sum
static inline void gk_abscissae(__float128 a, __float128 b, __float128* xgk) {

    const __float128* xgk_std = xgk_61;

    const __float128 center = 0.5Q * (b + a);
    const __float128 hl = 0.5Q * (b - a);

    for (int i = 0; i < DIV2(GK_POINTS-1); i += 2) {

        // Gauss points
        xgk[2*i] = center + hl * xgk_std[i+1];
        xgk[2*i+1] = center - hl * xgk_std[i+1];

        // Kronrod points
        xgk[2*i+2] = center + hl * xgk_std[i];
        xgk[2*i+3] = center - hl * xgk_std[i];

    }
    xgk[GK_POINTS-1] = center;

}

// Compute the weights (on a piece of the grid) for the 
// Gauss and Gauss-Kronrod methods
// NB: the same pattern is then used in gk_sum
static inline void gk_weigths(__float128 a, __float128 b, double* wgk) {

    const double* wgk_std = wgk_61;
    const double hl = 0.5 * (double)(b - a);

    for (int i = 0; i < DIV2(GK_POINTS-1); i += 2) {

        // Gauss-Kronrod points
        wgk[2*i]   = wgk_std[i+1] * hl;
        wgk[2*i+1] = wgk_std[i+1] * hl;
        wgk[2*i+2] = wgk_std[i]   * hl;
        wgk[2*i+3] = wgk_std[i]   * hl;

    }
    wgk[GK_POINTS-1] = wgk_std[DIV2(GK_POINTS-1)] * hl;

}

static inline void gk_weigths_gauss(__float128 a, __float128 b, double* wg) {

    const double* wg_std = wg_61;
    const double hl = 0.5 * (double)(b - a);

    for (int i = 0; i < DIV2(GK_POINTS-1); i += 2) {

        // Gauss points
        wg[2*i]   = wg_std[DIV2(i)] * hl;
        wg[2*i+1] = wg_std[DIV2(i)] * hl;
        wg[2*i+2] = 0.0;
        wg[2*i+3] = 0.0;

    }
    wg[GK_POINTS-1] = 0.0;

}

__float128* sl2cfoam_grid_uniform(int intervals) {

    __float128* grid;
    grid = (__float128*) malloc((intervals+1) * sizeof(__float128));

    __float128 dx = 1.0Q / intervals;

    grid[0] = 0.0Q;
    for (int i = 1; i < intervals; i++) {
        grid[i] = grid[i-1] + dx;
    }
    grid[intervals] = 1.0Q;

    return grid;

}

__float128* sl2cfoam_grid_harmonic(int intervals) {

    __float128* grid;
    grid = (__float128*) malloc((intervals+1) * sizeof(__float128));

    __float128 dx = 2.0Q / (intervals * (intervals+1));

    grid[0] = 0.0Q;
    for (int i = 1; i < intervals; i++) {
        grid[i] = grid[i-1] + i * dx;
    }
    grid[intervals] = 1.0Q;

    return grid;

}

__float128* sl2cfoam_gk_grid_abscissae(int intervals, __float128* grid) {

    __float128* xs;
    xs = (__float128*) malloc(GK_POINTS * intervals * sizeof(__float128));

    for (int i = 0; i < intervals; i++) {
        gk_abscissae(grid[i], grid[i+1], xs + GK_POINTS*i);
    }

    return xs;

}

double* sl2cfoam_gk_grid_weights(int intervals, __float128* grid) {

    double* ws;
    ws = (double*) malloc(GK_POINTS * intervals * sizeof(double));

    for (int i = 0; i < intervals; i++) {
        gk_weigths(grid[i], grid[i+1], ws + GK_POINTS*i);
    }

    return ws;

}

double* sl2cfoam_gk_grid_weights_gauss(int intervals, __float128* grid) {

    double* ws;
    ws = (double*) malloc(GK_POINTS * intervals * sizeof(double));

    for (int i = 0; i < intervals; i++) {
        gk_weigths_gauss(grid[i], grid[i+1], ws + GK_POINTS*i);
    }

    return ws;

}

double sl2cfoam_gk_grid(int intervals, double* ys, double* ms, 
                        double* wgks, double* wgs, double* abserr) {

    int nxs = GK_POINTS * intervals;
    
    long double res_kronrod = 0.0;

    // contract
    // TODO: verify if compensated summation (long double) is needed
    //       probably yes...
    for (int i = 0; i < nxs; i++) {
        res_kronrod += (long double)(ys[i] * ms[i] * wgks[i]);
    }

    if (wgs != NULL) {

        long double res_gauss = 0.0;
        for (int i = 0; i < nxs; i++) {
            res_gauss += (long double)(ys[i] * ms[i] * wgs[i]);
        }
        *abserr = (double)fabsl(res_kronrod-res_gauss);

    }

    return (double)res_kronrod;

}
