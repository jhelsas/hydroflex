#define hbarc (0.197) /* hbarc = 0.197 GeV fm*/
#define m_pi (0.139) /* m_pion = 0.139 GeV */
#define C_pi (0.0449427946273) /* (hbarc)*((45÷(3x128×π^2))^(1/3)) GeV fm */
#define C_qg (0.0194521040691) /* (hbarc)*((45÷(37x128×π^2))^(1/3)) GeV fm */ 
/*Checar se esse 37 esta certo, e se nao e 26.5 -- esse 37 esta certo, sao 7/8 mas tem que contar q e qbar
 * g_qg = 2x8 + (7/8)x 2(q qbar) x 2 (spin) x 3 (cor) x 2 (N_f) = 37
 * */
#define Sigma_pi (43.0308388315) /* (3xπ^2)/(90*(hbarc**3)) (GeV fm)^-3 */
#define Sigma_qg (530.713678922) /* (37xπ^2)/(90*(hbarc**3)) (GeV fm)^-3  */
/*#define EtaCharg_qg (0.778539493294)*/ /* (2/5)*(g_{3/2}(1)/g_{5/2}(1)) , ver Salinas 246-247 */
#define EtaCharg_qg (0.185104422545)  /* (15*zeta(3))/(pi^4) , ver wikipedia photon gas */

#define tau_R (1.0)
#define tau_R2 (0.5)
#define eta0 (1.0)

#define invtau_R (1.0)
#define inveta0tau_R2 (0.0)

int EoS_landau(SPHeq_particle *par);

int EoS_qg(SPHeq_particle *par);

int EoS_dust(SPHeq_particle *par);

int freezeout_const_time(int D,int N_sph,SPHeq_list *sph_eq,int Ny,int NpT,int Nphi,
                         int namingcase,char *hash,char *filepath);
