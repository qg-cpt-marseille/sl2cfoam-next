#include <quadmath.h>
#include <math.h>

#include "integration_qagp.h"
#include "utils.h"
#include "error.h"


// Norm function for errors in complex values.
static inline __float128 errnorm(__complex128 x) {
    return fabsq(crealq(x)) + fabsq(cimagq(x));
}

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
static const __float128 wg_61[15] =
{
    0.007968192496166605615465883474674Q,
    0.018466468311090959142302131912047Q,
    0.028784707883323369349719179611292Q,
    0.038799192569627049596801936446348Q,
    0.048402672830594052902938140422808Q,
    0.057493156217619066481721689402056Q,
    0.065974229882180495128128515115962Q,
    0.073755974737705206268243850022191Q,
    0.080755895229420215354694938460530Q,
    0.086899787201082979802387530715126Q,
    0.092122522237786128717632707087619Q,
    0.096368737174644259639468626351810Q,
    0.099593420586795267062780282103569Q,
    0.101762389748405504596428952168554Q,
    0.102852652893558840341285636705415Q
};

// weights for Kronrod points
static const __float128 wgk_61[31] =
{
    0.001389013698677007624551591226760Q,
    0.003890461127099884051267201844516Q,
    0.006630703915931292173319826369750Q,
    0.009273279659517763428441146892024Q,
    0.011823015253496341742232898853251Q,
    0.014369729507045804812451432443580Q,
    0.016920889189053272627572289420322Q,
    0.019414141193942381173408951050128Q,
    0.021828035821609192297167485738339Q,
    0.024191162078080601365686370725232Q,
    0.026509954882333101610601709335075Q,
    0.028754048765041292843978785354334Q,
    0.030907257562387762472884252943092Q,
    0.032981447057483726031814191016854Q,
    0.034979338028060024137499670731468Q,
    0.036882364651821229223911065617136Q,
    0.038678945624727592950348651532281Q,
    0.040374538951535959111995279752468Q,
    0.041969810215164246147147541285970Q,
    0.043452539701356069316831728117073Q,
    0.044814800133162663192355551616723Q,
    0.046059238271006988116271735559374Q,
    0.047185546569299153945261478181099Q,
    0.048185861757087129140779492298305Q,
    0.049055434555029778887528165367238Q,
    0.049795683427074206357811569379942Q,
    0.050405921402782346840893085653585Q,
    0.050881795898749606492297473049805Q,
    0.051221547849258772170656282604944Q,
    0.051426128537459025933862879215781Q,
    0.051494729429451567558340433647099Q
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

// Gaus-Kronrod summation given the function and measure
// evaluations.
static inline void gk_sum(__float128 a, __float128 b,
                          __complex128* ygk, __complex128* mgk,
                          __complex128* g_res, __complex128* gk_res) {

    const __float128* wg = wg_61;
    const __float128* wgk = wgk_61;

    const __float128 hl = 0.5Q * (b - a);

    __complex128 rg, rgk; 
    rg = 0.0Q; 
    rgk = ygk[GK_POINTS-1] * mgk[GK_POINTS-1] * wgk[DIV2(GK_POINTS-1)]; 

    int i;

    for (i = 0; i < DIV2(GK_POINTS-1); i += 2) {
        rg += (ygk[2*i]*mgk[2*i] + ygk[2*i+1]*mgk[2*i+1]) * wg[DIV2(i)]; // Gauss
    } 

    for (i = 0; i < DIV2(GK_POINTS-1); i += 2) { 
        rgk += (ygk[2*i]*mgk[2*i] + ygk[2*i+1]*mgk[2*i+1]) * wgk[i+1]; // Gauss
        rgk += (ygk[2*i+2]*mgk[2*i+2] + ygk[2*i+3]*mgk[2*i+3]) * wgk[i]; // Kronrod ext.
    } 

    *g_res = hl * rg; 
    *gk_res = hl * rgk; 

}

// macros for function evaluation
#define INTEGRATION_FEVAL(f, xs, ys, N) (*(f->function))(ys, xs, N, f->params)
#define INTEGRATION_MEVAL(f, xs, ys, N)  (*(f->measure))(ys, xs, N, f->params)

// Gauss-Kronrod 30-61 integration in interval (a b) with function evalution
static void qagp_gk(sl2cfoam_integration_function* f, __float128 a, __float128 b,
                    __complex128* gauss_res, __complex128* kronrod_res) {

    // compute abscissae
    __float128 xgk[GK_POINTS];
    gk_abscissae(a, b, xgk);

    // compute f at all points
    __complex128 ygk[GK_POINTS];
    INTEGRATION_FEVAL(f, xgk, ygk, GK_POINTS);

    // add measure function if defined
    __complex128 mgk[GK_POINTS];
    if (f->measure != NULL) {

        INTEGRATION_MEVAL(f, xgk, mgk, GK_POINTS);

    } else {

        // set measure to 1
        for (int i = 0; i < GK_POINTS; i++) {
            mgk[i] = 1.0Q;
        }

    }

    gk_sum(a, b, ygk, mgk, gauss_res, kronrod_res);

}

// contains the details of a subdivision step during
// adaptive integration
struct __subdiv {
    __float128 a;
    __float128 b;
    __complex128 ig;
    __complex128 ik;
    __float128 err;
    struct __subdiv* next;
    struct __subdiv* parent;
    int depth;
    long evals;
};

// The strategy of adaptive integration is the following:
// perform a first GK integration to get a rough estimate of the integral
// (which might be already good if function is simple and smooth)
// then call the main interval 'root' and start halving. At each halving
// estimate the error of both halves (difference btw gauss
// and kronrod estimates) and then move "left" or "right" to the interval
// with the largest estimated error. Assign the remaining evaluations
// proportionally (almost) to estimated error.
// At each step, if difference is small wrt global estimation,
// then accumulate the value and move "left" or "right", then do the same.
// Remaining evaluations are not lost but propagated.
// If there are no subintervals on the "side" then move up and try to move
// to the side again, and so on. At the end the last interval
// is reached, so stop there. Structures are freed in the traversing. 
sl2cfoam_qagp_retcode sl2cfoam_qagp(sl2cfoam_integration_function* f, double tol, __complex128* result, __float128* abserr,
                                    size_t* neval, __float128 dx_min, size_t max_eval) {

    size_t tot_eval = 0;
    const size_t eval_per_gk = 61;

    // get maximum allowed depth of subdivision
    const __float128* xgk_std = xgk_61;
    const int max_depth = max(0, (int)floorq(-log2q(dx_min / (1.0Q - xgk_std[0]))));
    if (max_depth < 4) {
        warning("max depth = %d too low, result may be inaccurate", max_depth)
    }

    // compute first approximation
    __complex128 ig, ik;
    qagp_gk(f, 0.0Q, 1.0Q, &ig, &ik);
    tot_eval += eval_per_gk;

    // estimation of global order of magnitude (Kronrod)
    __float128 is;
    is = errnorm(ik);

    struct __subdiv root;
    root.a = 0.0Q;
    root.b = 1.0Q;
    root.ig = ig;
    root.ik = ik;
    root.err = errnorm(ig - ik);
    root.next = NULL;
    root.parent = NULL;
    root.depth = 0;
    root.evals = max_eval;

    struct __subdiv* su = &root;
    struct __subdiv* su_next;
    struct __subdiv* su_up;
    struct __subdiv* sl;
    struct __subdiv* sr;

    __complex128 acc_int = 0.0Q;
    __float128 acc_abserr = 0.0Q;

    double fel, fer;
    int left_evals;
    bool max_depth_reached = false;
    for (;;) {

        check_interval:

        // check if error is within tolerance
        // or if maximum allowed depth has been reached
        // if yes, collect and go to next interval
        // if no, subdivide this interval
        if ((su->err < tol*is) || su->depth >= max_depth 
                               || su->evals < (int)eval_per_gk
                               || tot_eval > max_eval) {

            if (su->depth >= max_depth)
                max_depth_reached = true;

            // accumulate (kronrod)
            acc_int += su->ik;
            acc_abserr += su->err;

            // leftover evaluations to be passed over
            left_evals = imax(0, su->evals);
            
            if (su->next != NULL) {

                // move "left or right"
                su_next = su->next;
                free(su);
                su = su_next;
                su->evals += left_evals;

                goto check_interval;

            }

            while (su->parent != NULL) {

                // move "up"
                su_up = su->parent;
                free(su);
                su = su_up;

                if (su->next != NULL) {

                    // move "left or right"
                    su_next = su->next;
                    free(su);
                    su = su_next;
                    su->evals += left_evals;

                    goto check_interval;

                }

            }

            // reached last subinterval, stop here
            break;

        }

        // go on and subdivide this interval

        sl = (struct __subdiv*)malloc(sizeof(struct __subdiv));
        sr = (struct __subdiv*)malloc(sizeof(struct __subdiv));

        // subdivide
        sl->a = su->a;                         
        sl->b = sr->a = (su->a + su->b) * 0.5Q;
        sr->b = su->b;                         
                   
        sl->parent = su;                       
        sr->parent = su;                       
        sl->depth = sr->depth = su->depth + 1;

        // compute integrals and errors
        qagp_gk(f, sl->a, sl->b, &(sl->ig), &(sl->ik));
        qagp_gk(f, sr->a, sr->b, &(sr->ig), &(sr->ik));
        sl->err = errnorm(sl->ig - sl->ik);
        sr->err = errnorm(sr->ig - sr->ik);

        su->evals -= 2 * eval_per_gk;
        tot_eval += 2 * eval_per_gk;

        // assign remaining evaluations proportionally
        // to estimated error
        // use f(x) = sqrt(x) to "give less importance" to error
        // (since it is an estimate)
        fel = sqrt((double)sl->err);
        fer = sqrt((double)sr->err);
        sl->evals = (long)(su->evals * fel / (fel + fer));
        sr->evals = (long)(su->evals * fer / (fel + fer));

        // move "left or right"
        if (sl->err >= sr->err) {
            sl->next = sr;                         
            sr->next = NULL;    
            su = sl;
        } else {
            sl->next = NULL;                         
            sr->next = sl;    
            su = sr;
        }

    }

    // set values computed till here
    *result = acc_int;
    *abserr = acc_abserr;
    *neval = tot_eval;

    sl2cfoam_qagp_retcode ret;
    ret = QAGP_SUCCESS;

    // add error flags if there were errors

    if (tot_eval >= max_eval)
        ret |= QAGP_MAX_EVAL;

    if (max_depth_reached)
        ret |= QAGP_MAX_DEPTH;

    return ret;

}
