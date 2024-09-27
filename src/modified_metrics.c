#include "modified_metrics.h"

#include "decs.h"
#include "geometry.h"
#include "debug_tools.h"
//personal note, try to implement max EH 
int theory = 1; //Theory of Gravity: 0. for General Relativity, 1. for EdGB, and 2. for DCS
double zeta = 0.1; //Deviation from Gravity: Should be between 0. and 0.3
int cano_solution = 0; //choose 1 to use the cano solution for the modified metrics or 0 for the regular solution
//the main difference is that the cano solution is up to 31st order in spin and therefore takes much longer to run
//the regular solution is only up to 5th order in spin

/*Because of the length of the cano solutions it's likely easiest to navigate this file via searching for the function you want*/

double a;
double Rh;

void gcov_EdGB_ks(double r, double th, double gcov[NDIM][NDIM]){
//Return the Einstein dilaton Gauss Bonet metric in KS coordinates
//In order to make plugging in new metrics easier I have changed the typical shorthand variables we use
//the suffices have the following meanings: -sq=squared, -cub=cubed, -hcub= hypercube (^4), -surs=sursolid (^5)
  double u=cos(th);
  double u2 = u*u;
  double u4 = u2*u2;
  double u6 = u2*u2*u2;
  double M[NDIM][NDIM];
  double gmod[NDIM][NDIM];
  double expM[NDIM][NDIM];
  double expMT[NDIM][NDIM];
  double MT[NDIM][NDIM];
  double f;

  MUNULOOP{
    gcov[mu][nu] = 0.;
    M[mu][nu]=0.;
    gmod[mu][nu]=0.;
    expM[mu][nu]=0.;
    expMT[mu][nu]=0.;
    MT[mu][nu]=0.;
  } 

  if(r<=Rh){
    fprintf(stderr,"trying to go inside the EH!!! \n");
    M[0][0]= r*(0.09881663310865521 - 0.06589959486016632*u2 + 
      0.09483604081290098*u4 + r*(-0.0508082639442996 + 
        0.03603450867902285*u2 - 0.04755871830571102*u4));
    M[0][1]= r*(1.1660157755065736 - 0.923446854438797*u2 + 0.4445882546481567*u4 + 
      r*(-0.5358018673193612 + 0.45031722475873553*u2 - 
        0.22156152685943775*u4));
    M[0][2]= 0.;
    M[0][3]=r*(0.26629173567914655 - 0.012224165629961845*u2 - 
      0.4504719467791443*u4 + 0.19640437672995958*u6 + 
      r*(-0.117216423840705 - 0.007938446653934146*u2 + 
        0.22304710326845398*u4 - 0.09789223277381481*u6));

    M[1][0]=0.;
    M[1][1]=r*(-0.21770433903942524 + 0.18701808101289774*u2 - 
      0.09675918672289041*u4 + r*(0.1031193586421407 - 
        0.09134547021704825*u2 + 0.04792976859542179*u4));
    M[1][2]=0.;
    M[1][3]=r*(-0.05211790724575234 + 0.10282977052005256*u2 - 
      0.04878574904287365*u4 - 0.0019261142314265672*u6 + 
      r*(0.024949669974718935 - 0.04977020053092391*u2 + 
        0.024521350206178745*u4 + 0.0002991803500262339*u6));

    M[2][0]=0.;
    M[2][1]=0.;
    M[2][2]=(r*(-0.2456694322274893 + 0.11802793355278088*u2 - 
      0.06388746774868627*u4 + r*(0.11088665306540471 - 
        0.06020804162310437*u2 + 0.033009174906014245*u4)))*sin(th)*sin(th);
    M[2][3]=0.;

    M[3][0]=0.;
    M[3][1]=0.;
    M[3][2]=0.;
    M[3][3]=r*(-0.10427504873012293 - 0.05238594045593648*u2 - 
      0.034867977237335264*u4 + r*(0.04443039634840591 + 
        0.02057044808324951*u2 + 0.01868694191665917*u4));

    f = r*(1.0629450418611497 + 0.01171271935304809*u2 - 
    0.00493864955500953*u4 + r*(-0.01745401695228727 - 
      0.0032478173991280367*u2 + 0.0013694370597876334*u4));
  }else if(r>Rh){
        M[0][0]=(4.291164274322169*u4 + pow(r,3.)*(-0.2109264049823971 - 
        0.727527995711403*u2 - 0.2289082277890738*u4) + 
      pow(r,2.)*(0.003213989800925891 - 2.8729270695242883*u2 - 
        0.07200152505487534*u4) + pow(r,4.)*(0.7091703981563269 + 
        0.03414119433807664*u2 - 0.039940983006676915*u4) + 
      pow(r,8.)*(-0.056750235572314435 - 0.010559999053525839*u2 + 
        0.00445260686733901*u4) + pow(r,7.)*(-0.1071071860501588 + 
        0.02587578792636418*u2 + 0.00890521373467802*u4) + 
      pow(r,5.)*(0.3036965654525125 + 0.35876703857586284*u2 + 
        0.016769533525680457*u4) + pow(r,6.)*(0.16405391299415248 + 
        0.14311498384751617*u2 + 0.02155839824857682*u4) + 
      r*(0.25435376653421765*u2 + 0.8575513355876029*u4))/pow(r,10.);
        M[0][1]=(6.788644338118022*u4 + pow(r,8.)*(0.1847953257927588 - 0.021119998107051677*
         u2 + 0.00890521373467802*u4) + pow(r,7.)*(0.5012698809657435 - 
        0.017972840146824282*u2 + 0.01781042746935604*u4) + 
      pow(r,6.)*(1.1678820823398104 - 0.08665926837000991*u2 + 
        0.04311679649715364*u4) + pow(r,5.)*(1.3968855321134972 - 
        0.5214758738558606*u2 + 0.05509367847082749*u4) + 
      pow(r,4.)*(1.2803751437389894 - 1.9801727988754083*u2 + 
        0.11100987467318106*u4) + pow(r,3.)*(-0.19674673106180973 - 
        3.631607281837495*u2 + 0.4584931004015606*u4) + 
      pow(r,2.)*(0.027577060783827698 - 4.600749569082939*u2 + 
        1.927572690712894*u4) + r*(0.30636166840678114*u2 + 
        4.428944923865755*u4))/pow(r,10.);
        M[0][2]=0.;
        M[0][3]=(-3.394322169059011*u4 + 3.394322169059011*u6 + 
      pow(r,8.)*(0.12158613671972822 - 0.11102613766620238*u2 - 
        0.01501260592086485*u4 + 0.00445260686733901*u6) + 
      pow(r,7.)*(0.26003646165384053 - 0.27131451731494105*u2 + 
        0.0023728419264224393*u4 + 0.00890521373467802*u6) + 
      pow(r,6.)*(0.15034750415271614 - 0.18684369578994947*u2 + 
        0.014937793388656521*u4 + 0.02155839824857682*u6) + 
      pow(r,5.)*(0.13078293898666954 - 0.09663728782559695*u2 - 
        0.05834564172670751*u4 + 0.024199990565634915*u6) + 
      pow(r,4.)*(-0.32370570453399117 + 0.8745571659528123*u2 - 
        0.5881531198749078*u4 + 0.0373016584560868*u6) + 
      pow(r,3.)*(0.3444360776914127 + 0.8053120998423906*u2 - 
        1.3136486338971856*u4 + 0.16390045636338207*u6) + 
      pow(r,2.)*(0.008855468153400483 + 2.0874440102442646*u2 - 
        2.8957358880860653*u4 + 0.7994364096883999*u6) + 
      r*(-0.3627447439778267*u2 - 1.5398496895507068*u4 + 
        1.9025944335285334*u6))/pow(r,10.);

        M[1][0] = 0.;
        M[1][1] = (-2.497480063795853*u4 + r*(-0.052007901872563524*u2 - 
        0.7819199040676263*u4) + pow(r,2.)*(-0.02436307098290181 + 
        1.7686911221641155*u2 - 0.17276456721199093*u4) + 
      pow(r,5.)*(-0.33186954053463574 - 0.04713633522918607*u2 - 
        0.028384689943681703*u4) + pow(r,6.)*(-0.22196976327416232 - 
        0.025450775001065937*u2 - 0.02446218735307713*u4) + 
      pow(r,4.)*(-0.5497817903708182 + 0.27738723097161644*u2 - 
        0.010676666072864816*u4) + pow(r,7.)*(-0.017878650978800745 + 
        0.01361897560289736*u2 - 0.00890521373467802*u4) + 
      pow(r,8.)*(0.0018711207888258974 + 0.022693577087165375*u2 - 
        0.00445260686733901*u4) + pow(r,3.)*(-0.011217539355069367 + 
        0.7538723496096127*u2 + 0.04135389749183582*u4))/pow(r,10.);
        M[1][2]=0.;
        M[1][3]=(0.10478195488721805*u2 - 0.2095639097744361*u4 + 
      0.10478195488721805*u6 + pow(r,4.)*(-0.07454768676004499 + 
        0.14076463686934715*u2 - 0.05788621345855937*u4 - 
        0.008330736650742812*u6) + pow(r,5.)*(-0.03306014626262768 + 
        0.05858256643029559*u2 - 0.017984694072708127*u4 - 
        0.007537726094959777*u6) + pow(r,6.)*(-0.014357782202486562 + 
        0.023325964842828922*u2 - 0.003890913976466552*u4 - 
        0.005077268663875806*u6) + pow(r,3.)*(-0.10234633569089166 + 
        0.2020710498436015*u2 - 0.09710309261452797*u4 - 
        0.0026216215381818465*u6) + pow(r,7.)*(-0.0033719892357476668 + 
        0.004158778725804516*u2 + 0.0014395139436126551*u4 - 
        0.002226303433669505*u6) + pow(r,8.)*(0.004042832787549678 - 
        0.009322832314312596*u2 + 0.007506302960432425*u4 - 
        0.002226303433669505*u6) + pow(r,2.)*(-0.09979646693535539 + 
        0.22353023569085523*u2 - 0.14767107057564424*u4 + 
        0.02393730182014443*u6) + r*(0.0015552425331560671 + 
        0.06452527650736407*u2 - 0.13371628061419635*u4 + 
        0.0676357615736762*u6))/pow(r,9.);

        M[2][0] = 0.;
        M[2][1] = 0.;
        M[2][2] = ((-2.497480063795853*u4 + pow(r,8.)*(-0.00041852575015099877 + 
        0.012133578033639538*u2) + r*(-0.052007901872563524*u2 - 
        0.7819199040676263*u4) + pow(r,2.)*(0.0013914126287246588 + 
        1.7171821549408626*u2 - 0.14701008360036444*u4) + 
      pow(r,6.)*(-0.20399052197974238 - 0.061409257589905675*u2 - 
        0.006482946058657265*u4) + pow(r,5.)*(-0.30359291081514894 - 
        0.10368959466815979*u2 - 0.00010806022419483593*u4) + 
      pow(r,7.)*(-0.0072303365947720705 - 0.007939422645679456*u2 + 
        0.001249323593073593*u4) + pow(r,9.)*(-0.056750235572314435 - 
        0.010559999053525839*u2 + 0.00445260686733901*u4) + 
      pow(r,4.)*(-0.5106316084099219 + 0.19908686704982395*u2 + 
        0.028473515888031427*u4) + pow(r,3.)*(0.028197926131824688 + 
        0.6750414186358245*u2 + 0.08076936297872989*u4))/pow(r,10.))*sin(th)*sin(th);
        M[2][3]= 0.;

        M[3][0] = 0.;
        M[3][1] =0.;
        M[3][2] =0.;
        M[3][3] =(0.011715052283488538*pow(r,8.) - 2.497480063795853*u4 + 
      r*(-0.2615718116469996*u2 - 0.5723559942931902*u4) + 
      pow(r,6.)*(-0.11989098816006717 - 0.1484125805140812*u2 - 
        0.003579156954156954*u4) + pow(r,7.)*(0.032264426934489464 - 
        0.04743418617494099*u2 + 0.001249323593073593*u4) + 
      pow(r,9.)*(-0.056750235572314435 - 0.010559999053525839*u2 + 
        0.00445260686733901*u4) + pow(r,5.)*(-0.12622090757557206 - 
        0.29267675432573786*u2 + 0.011507096193806405*u4) + 
      pow(r,2.)*(0.02403541117403899 + 1.5335121496365696*u2 + 
        0.014015923158614446*u4) + pow(r,4.)*(-0.2667887550672424 - 
        0.07866292517738821*u2 + 0.06238045477256399*u4) + 
      pow(r,3.)*(0.2672063254894295 + 0.3487429501510368*u2 + 
        0.1680594321059128*u4))/pow(r,10.);

        f= 0.056750235572314435 + r + 0.01055999905352584*u2 - 0.00445260686733901*u4;
        if(r>1.0E15 || r<-1.0E15){
          fprintf(stderr,"big r: %.15lf \n corresponding f: %.15lf \n",r,f);
        }
    }

   // double cth = cos(th);
   // double sth = sin(th);
  //double ssq = sth * sth;
    double rhosq = f * f + a * a * u2;
    double sth = sin(th);
    double s2 = sth*sth;
    
    double gcov_kerr[NDIM][NDIM];
    gcov_kerr[0][0]=-1. + (2.*f)/rhosq;
    gcov_kerr[0][1]=(2.*f)/rhosq;
    gcov_kerr[0][2]=0.0;
    gcov_kerr[0][3]=-2.*f*a*s2/rhosq;
    

    gcov_kerr[1][0]=gcov_kerr[0][1];
    gcov_kerr[1][1]=1.+ (2.*f)/rhosq;
    gcov_kerr[1][2]=0.0;
    gcov_kerr[1][3]=(-a*s2*(1.+(2.*f)/rhosq));

    gcov_kerr[2][0]=0.0;
    gcov_kerr[2][1]=0.0;
    gcov_kerr[2][2]=rhosq;
    gcov_kerr[2][3]=0.0;

    gcov_kerr[3][0]=gcov_kerr[0][3];
    gcov_kerr[3][1]=gcov_kerr[1][3];
    gcov_kerr[3][2]=gcov_kerr[2][3];
    gcov_kerr[3][3] =s2*(rhosq+ a*a*s2*(1.+(2.*f)/rhosq));
   
    
    
    MUNULOOP{
      MT[nu][mu] = M[mu][nu];
    }
    matrix_exponential(MT,expMT);
    matrix_multiply(expMT,gcov_kerr,gcov);

    matrix_exponential(M,expM);
    matrix_multiply(gcov,expM,gcov);
    //fprintf(stderr,"r: %.15lf,th: %.15lf \n",r,th);
    //fprintf(stderr,"%.15lf,%.15lf,%.15lf%.15lf \n %.15lf,%.15lf,%.15lf,%.15lf \n %.15lf,%.15lf,%.15lf,%.15lf \n %.15lf,%.15lf,%.15lf,%.15lf \n",gcov[0][0],gcov[0][1],gcov[0][2],gcov[0][3],gcov[1][0],gcov[1][1],gcov[1][2],gcov[1][3],gcov[2][0],gcov[2][1],gcov[2][2],gcov[2][3],gcov[3][0],gcov[3][1],gcov[3][2],gcov[3][3]);
    
    //print_matrix("gcov EdGB",gcov);
    //print_matrix("gcov GR",gcov_kerr);
  
}


void gcov_DCS_ks(double r, double th, double gcov[NDIM][NDIM]){
//Return the Dynamical Chern-Simons metric in KS coordinates
  double u=cos(th);
  double u2 = u*u;
  double u4 = u2*u2;
  double u6 = u2*u2*u2;
  double M[NDIM][NDIM];
  double gmod[NDIM][NDIM];
  double expM[NDIM][NDIM];
  double MT[NDIM][NDIM];
  double f;
  double expMT[NDIM][NDIM];


  MUNULOOP{
    gcov[mu][nu] = 0.;
    M[mu][nu]=0.;
    gmod[mu][nu]=0.;
    expM[mu][nu]=0.;
    expMT[mu][nu]=0.;
    MT[mu][nu]=0.;
  } 

  if(r<=Rh){
       fprintf(stderr,"try to go inside EH!\n");
        M[0][0]=r*(-0.10186551426782008 + 0.2044681509464842*u2 - 0.047577127121860865*u4 + r*(0.045274692542072226 - 
        0.09476108773453298*u2 + 0.024336657157255144*u4));
        M[0][1]= r*(0.045412209337770834 + 0.6948770087146782*u2 - 0.22740730325285005*u4 + r*(-0.018243440528231694 - 
        0.3195382642179189*u2 + 0.11136116213415795*u4));
        M[0][2]= 0.;
        M[0][3]= r*(0.9245158222592134 - 1.5540867827958567*u2 + 0.7817181596052425*u4 - 0.15214719906859933*u6 + r*(-0.39831379601518213 + 
        0.6873438631794827*u2 - 0.3629692425201996*u4 + 
        0.07393917535589899*u6));


        M[1][0] = 0.;
        M[1][1] = r*(0.016459528849326405 - 0.0877649193211191*u2 + 0.013672933340417168*u4 + r*(-0.005981133205508223 + 
        0.0397302074948605*u2 - 0.008047529144189837*u4));
        M[1][2]=0.;
        M[1][3]=r*(0.08767766008723328 - 0.03217485989886997*u2 - 0.04878094127510386*u4 - 0.006721858913259437*u6 + r*(-0.027339833822029428 + 0.0032782178868501132*u2 + 
        0.021849341728100158*u4 + 0.002212274207079146*u6));

        M[2][0] = 0.;
        M[2][1] =0. ;
        M[2][2] =(r*(-0.005891088388619743 - 0.08213157970162656*u2 + 0.04869729617963716*u4 + r*(0.0012258296731578948 + 
        0.04059152589803598*u2 - 0.022648614125266385*u4)))*sin(th)*sin(th);
        M[2][3]=0.;

        M[3][0] = 0.;
        M[3][1] = 0.;
        M[3][2] = 0.;
        M[3][3] = r*(0.0779005273559132 - 0.18750384045886062*u2 + 0.07027794119233827*u4 + r*(-0.03758218403345795 + 
        0.08957607617157008*u2 - 0.03282515069218463*u4));
        
        //f is the radius correction we need to pass to kerr metric
        f =r*(1.020691229419899 - 0.031030208008011838*u2 - 0.01092031131121305*u4 + r*(-0.005485265354929308 + 0.008226138789940869*u2 + 0.0028949853140593047*u4));
  }else if(r>Rh){
        M[0][0]=(-3.764802631578948*u4 + r*(0.20737746740621577*u2 - 
        1.2411199826309292*u4) + pow(r,2.)*(0.004397491677732722 + 
        2.275428565488274*u2 - 0.25329799985773566*u4) + 
        pow(r,8.)*(-0.019512591789299767 + 0.02926262957652467*u2 + 
        0.010298255966503736*u4) + pow(r,6.)*(-0.07871641654069922 + 
        0.018403555854965288*u2 + 0.011519474758872082*u4) + 
        pow(r,5.)*(-0.1601638195572365 + 0.011472388631105742*u2 + 
        0.01770559700437244*u4) + pow(r,7.)*(-0.03924720003976625 + 
        0.02044660472170251*u2 + 0.02059651193300747*u4) + 
        pow(r,4.)*(-0.2062015656196536 + 0.3715397076125874*u2 + 
        0.030607798809265254*u4) + pow(r,3.)*(-0.19117702886935253 + 
        1.0350857068689334*u2 + 0.07147163859562965*u4))/pow(r,10.);
        
        M[0][1]=(-5.950657894736843*u4 + r*(0.20903884720779903*u2 - 
        4.565478341021519*u4) + pow(r,2.)*(0.030794723900581383 + 
        3.383010068231072*u2 - 2.1399936508080724*u4) + 
        pow(r,3.)*(-0.1496858869925613 + 3.1029202972978545*u2 - 
        0.5753275196143051*u4) + pow(r,4.)*(-0.09506265415433443 + 
        1.8595028809086818*u2 - 0.07821619255002231*u4) + 
        pow(r,8.)*(0.03680861200741224 + 0.05852525915304934*u2 + 
        0.02059651193300747*u4) + pow(r,5.)*(-0.012651581837955407 + 
        0.7163988910916779*u2 + 0.02211495126091475*u4) + 
        pow(r,6.)*(0.07232253531051423 + 0.22249925201839119*u2 + 
        0.023038949517744163*u4) + pow(r,7.)*(0.05973819193995381 + 
        0.06827691730484674*u2 + 0.04119302386601494*u4))/pow(r,10.);
        M[0][2]=0.;
        M[0][3]= (2.9753289473684217*u4 - 2.9753289473684217*u6 + 
        r*(-0.3301529518690068*u2 + 3.1596148319314237*u4 - 
        2.8294618800624165*u6) + pow(r,2.)*(0.008578113558176309 - 
        1.9548229788672007*u2 + 3.4050066693919323*u4 - 
        1.458761804082908*u6) + pow(r,3.)*(0.3315641345145541 - 
        2.9057414924554914*u2 + 3.0361396822695923*u4 - 
        0.461962324328655*u6) + pow(r,4.)*(0.39370132858720774 - 
        2.167692591796219*u2 + 1.857264696904211*u4 - 
        0.08327343369519986*u6) + pow(r,5.)*(0.9987748278913282 - 
        1.840302530914861*u2 + 0.8372567620733599*u4 + 
        0.004270940950172779*u6) + pow(r,8.)*(0.3226400239145535 - 
        0.3519026534910782*u2 + 0.018964373610020936*u4 + 
        0.010298255966503736*u6) + pow(r,6.)*(0.9441539653309928 - 
        1.222496704008794*u2 + 0.2668232639189292*u4 + 
        0.011519474758872082*u6) + pow(r,7.)*(0.7174242038941537 - 
        0.7898827927584792*u2 + 0.05186207693131773*u4 + 
        0.02059651193300747*u6))/pow(r,10.);

        M[1][0] = 0.;
        M[1][1] =  (2.185855263157895*u4 + pow(r,3.)*(-0.03416402565299624 - 
        0.48939094433013075*u2 - 0.0694329299472212*u4) + 
        pow(r,4.)*(-0.029868113792055086 - 0.12779951461553218*u2 - 
        0.049418835957736065*u4) + pow(r,5.)*(-0.015665261737899675 + 
        0.07028917056179607*u2 - 0.029196090158295766*u4) + 
        pow(r,7.)*(0.006093802057191114 + 0.0060208618032693405*u2 - 
        0.02059651193300747*u4) + pow(r,6.)*(-0.008285434801374236 + 
        0.038395483503946296*u2 - 0.014392098047352914*u4) + 
        pow(r,8.)*(0.040556426561124084 - 0.05364943007715063*u2 - 
        0.010298255966503736*u4) + pow(r,2.)*(-0.026397232222848656 - 
        1.1827213707401538*u2 + 0.16206700830456913*u4) + 
        r*(-0.0016613798015832518*u2 + 0.7709044110221688*u4))/pow(r,10.);
        M[1][2]= 0.;
        M[1][3]=(0.11281676413255362*u2 - 0.22563352826510724*u4 + 
        0.11281676413255362*u6 + pow(r,4.)*(-0.07720323170405753 + 
        0.14576209104746674*u2 - 0.05991448698276085*u4 - 
        0.008644372360648345*u6) + pow(r,5.)*(-0.03212746825771086 + 
        0.056741178976613044*u2 - 0.017099953180093504*u4 - 
        0.0075137575388086785*u6) + pow(r,6.)*(-0.010353864572767907 + 
        0.017827860455817794*u2 - 0.002121331434403294*u4 - 
        0.005352664448646592*u6) + pow(r,8.)*(0.141807420167977 - 
        0.1271761053797147*u2 - 0.009482186805010468*u4 - 
        0.005149127983251868*u6) + pow(r,7.)*(0.013633334595171994 - 
        0.011195420057222646*u2 + 0.002711213445302516*u4 - 
        0.005149127983251868*u6) + pow(r,3.)*(-0.10714813519542485 + 
        0.21170106044976475*u2 - 0.10195771531325491*u4 - 
        0.002595209941084969*u6) + pow(r,2.)*(-0.10346977678554554 + 
        0.234990701941892*u2 - 0.15957207352714736*u4 + 
        0.02805114837080091*u6) + r*(0.0012855403069909927 + 
        0.07377655505803313*u2 - 0.15140973103703922*u4 + 
        0.07634763567201511*u6))/pow(r,9.);

        M[2][0] = 0.;
        M[2][1] = 0.;
        M[2][2] = ((2.185855263157895*u4 + pow(r,8.)*(0.0033429657288016426 - 
        0.024386800500625968*u2) + pow(r,3.)*(0.006730793813251797 - 
        0.5711805832626268*u2 - 0.028538110480973156*u4) + 
        pow(r,7.)*(0.006907947918770011 - 0.005498612955602743*u2 - 
        0.009891183035714287*u4) + pow(r,4.)*(0.010639314047962048 - 
        0.2088143702955665*u2 - 0.008911408117718933*u4) + 
        pow(r,5.)*(0.013113976137320347 + 0.012730694811356026*u2 - 
        0.00041685228307574617*u4) + pow(r,6.)*(0.009614703564723954 + 
        0.0025952067717499227*u2 + 0.0035080403187452754*u4) + 
        pow(r,9.)*(-0.019512591789299767 + 0.02926262957652467*u2 + 
        0.010298255966503736*u4) + pow(r,2.)*(0.000149323899600326 - 
        1.2358144829850521*u2 + 0.18861356442701813*u4) + 
        r*(-0.0016613798015832518*u2 + 0.7709044110221688*u4))/pow(r,10.))*sin(th)*sin(th);
        M[2][3]= 0.;

        M[3][0] = 0.;
        M[3][1] = 0.;
        M[3][2] = 0.;
        M[3][3] = (-0.021043834771824328*pow(r,8.) + 2.185855263157895*u4 + 
        pow(r,7.)*(0.03337541444374187 - 0.031966079480574595*u2 - 
        0.009891183035714287*u4) + pow(r,6.)*(0.09176977844624386 - 
        0.08243249139825083*u2 + 0.006380663607226107*u4) + 
        pow(r,9.)*(-0.019512591789299767 + 0.02926262957652467*u2 + 
        0.010298255966503736*u4) + pow(r,5.)*(0.1962996774206555 - 
        0.18194549962590242*u2 + 0.011073640870847584*u4) + 
        pow(r,4.)*(0.2654430122788289 - 0.4989350764842804*u2 + 
        0.026405599840128257*u4) + pow(r,3.)*(0.2545651668505909 - 
        0.9160120725078158*u2 + 0.0684590057268767*u4) + 
        pow(r,2.)*(0.024124799408067325 - 1.4390317859599984*u2 + 
        0.36785539189349736*u4) + r*(-0.2272949080666905*u2 + 
        0.9965379392872759*u4))/pow(r,10.);

        f= 0.01951259178929976 + r - 0.02926262957652467*u2 - 0.010298255966503736*u4;
        
  }

   // double cth = cos(th);
   // double sth = sin(th);
  //double ssq = sth * sth;
    double rhosq = f * f + a * a * u2;
    double sth = sin(th);
    double s2 = sth*sth;
    
    
    double gcov_kerr[NDIM][NDIM];
    gcov_kerr[0][0]=-1. + (2.*f)/rhosq;
    gcov_kerr[0][1]=(2.*f)/rhosq;
    gcov_kerr[0][2]=0.0;
    gcov_kerr[0][3]=(-2.*f*a*s2)/rhosq;

    gcov_kerr[1][0]=gcov_kerr[0][1];
    gcov_kerr[1][1]=1.+ (2.*f)/rhosq;
    gcov_kerr[1][2]=0.0;
    gcov_kerr[1][3]=(-a*s2*(1.+(2.*f)/rhosq));

    gcov_kerr[2][0]=0.0;
    gcov_kerr[2][1]=0.0;
    gcov_kerr[2][2]=rhosq;
    gcov_kerr[2][3]=0.0;

    gcov_kerr[3][0]=gcov_kerr[0][3];
    gcov_kerr[3][1]=gcov_kerr[1][3];
    gcov_kerr[3][2]=gcov_kerr[2][3];
    gcov_kerr[3][3] =s2*(rhosq+ a*a*s2*(1.+(2.*f)/rhosq));
   
    
    
    MUNULOOP{
      MT[nu][mu] = M[mu][nu];
    }
    matrix_exponential(MT,expMT);
    matrix_multiply(expMT,gcov_kerr,gcov);

    matrix_exponential(M,expM);
    matrix_multiply(gcov,expM,gcov);
   // fprintf(stderr,"r: %.15lf,th: %.15lf \n",r,th);
    
   // print_matrix("gcov dCS",gcov);
   // print_matrix("gcov GR",gcov_kerr);

}


double get_EdGB_Event_Horizon(){
    //returns the largest event horizon for EdGB 
    //solely leaving this function existing for whenever we do the general case and have an actual expression to use
    return 1.8031677760259373;
}


double get_dCS_Event_Horizon(){
    //returns the largest event horizon for dCS for theta between 0. and pi
    //solely leaving this function existing for whenever we do the general case and have an actual expression to use
    return 1.8860736975381673;
}

double event_horizon(){
  //return Rh;
  if(theory==0){
    return 1. + sqrt(1. - a * a);
  }else if(theory==1){
    //return 1. + sqrt(1. - a * a);
    return  1.8031677760259373;
  }else if(theory==2){
    return 1.8860736975381673;
    //return get_dCS_Event_Horizon();
  }
}


// Function to multiply two matrices
void matrix_multiply(double A[NDIM][NDIM], double B[NDIM][NDIM], double result[NDIM][NDIM]) {
    double temp[NDIM][NDIM];
    MUNULOOP temp[mu][nu]=0.;
    
    MUNULOOP{
      for (int k = 0; k < NDIM; k++) {
                temp[mu][nu] += A[mu][k] * B[k][nu];
      }
    }
    // Copy the result to the output matrix
    MUNULOOP result[mu][nu] = temp[mu][nu];
}



// Function to calculate the matrix exponential using 4th-order Taylor approximation
void matrix_exponential(double A[NDIM][NDIM], double expA[NDIM][NDIM]) {
    
    MULOOP expA[mu][mu]=1.0;

    // Temporary matrices for powers of A
    double A_pow[NDIM][NDIM] = {0};
    double temp[NDIM][NDIM] = {0};
    
    // expA = I + A
    MUNULOOP expA[mu][nu]=A[mu][nu]+expA[mu][nu];

    // A^2 / 2!
    matrix_multiply(A, A, A_pow);
    MUNULOOP temp[mu][nu] = A_pow[mu][nu]/2.;
    MUNULOOP expA[mu][nu]=temp[mu][nu]+expA[mu][nu];

    // A^3 / 3!
    matrix_multiply(A_pow, A, A_pow);
    MUNULOOP temp[mu][nu] = A_pow[mu][nu]/6.;
    MUNULOOP expA[mu][nu]=temp[mu][nu]+expA[mu][nu];

    // A^4 / 4!
    matrix_multiply(A_pow, A, A_pow);
    MUNULOOP temp[mu][nu] = A_pow[mu][nu]/24.;
    MUNULOOP expA[mu][nu]=temp[mu][nu]+expA[mu][nu];
}

