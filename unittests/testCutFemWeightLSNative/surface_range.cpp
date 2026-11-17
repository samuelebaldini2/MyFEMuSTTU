// Standalone assembled Q2 velocity residual. Uses the production moment fitter,
// FE basis and Gauss rules; this is NOT the PETSc/FEMuS system residual.
#include "CutFemWeight_LS.hpp"
#include "TangentialSurfaceIntegration.hpp"

#include "Basis.hpp"
#include <fstream>
#include <iomanip>
#include <map>
#include <array>
#include <chrono>
#include <sstream>
#include <memory>
using Point=std::array<double,3>;
struct Cell { unsigned id,child=0; std::vector<Point> x; };
femus::HexBiquadratic hex;
double sphere(Point x) { double g=.04; for(double a:x) g-=(a-.5)*(a-.5); return g; }
double phi(unsigned i,const Point& p) { return hex.eval_phi(hex.GetIND(i),p.data()); }
Point grad(unsigned i,const Point& p) { return {hex.eval_dphidx(hex.GetIND(i),p.data()),hex.eval_dphidy(hex.GetIND(i),p.data()),hex.eval_dphidz(hex.GetIND(i),p.data())}; }
Point map(const Cell& c,Point p) { Point x={}; for(unsigned i=0;i<27;++i) for(unsigned d=0;d<3;++d) x[d]+=phi(i,p)*c.x[i][d]; return x; }
std::vector<Cell> mesh(std::string geometry,std::string path) {
 std::vector<Cell> cells;
 if(geometry=="B") {
  std::ifstream in(path); if(!in) throw std::runtime_error("Mesh not found");
  std::string s; while(in>>s && s!="NDFVL"){} unsigned nn,ne,tmp; in>>nn>>ne>>tmp>>tmp>>tmp>>tmp;
  while(in>>s && s!="COORDINATES"){} in>>s; std::vector<Point> nodes(nn+1);
  for(unsigned i=0;i<nn;++i){unsigned id;in>>id;for(double& a:nodes[id]) in>>a;}
  while(in>>s && s!="ELEMENTS/CELLS"){} in>>s;
  const unsigned permutation[27]={4,16,0,15,23,11,7,19,3,12,20,8,25,26,24,14,22,10,5,17,1,13,21,9,6,18,2};
  for(unsigned e=0;e<ne;++e){Cell c;unsigned type,n;in>>c.id>>type>>n;if(n!=27)throw std::runtime_error("Not HEX27"); c.x.resize(27);for(unsigned i=0;i<27;++i){unsigned j;in>>j;c.x[permutation[i]]=nodes.at(j);} if(std::getenv("CUTFEM_PROJECT_Q1")){auto corners=c.x;for(unsigned n=0;n<27;++n){auto p=hex.GetXcoarse(n);Point x={};for(unsigned v=0;v<8;++v){auto t=hex.GetXcoarse(v);double w=(1+t[0]*p[0])*(1+t[1]*p[1])*(1+t[2]*p[2])/8;for(unsigned d=0;d<3;++d)x[d]+=w*corners[v][d];}c.x[n]=x;}} cells.push_back(c);}
 } else {
  for(unsigned k=0;k<4;++k) for(unsigned j=0;j<4;++j) for(unsigned i=0;i<4;++i) {
   Cell c; c.id=cells.size()+1;c.x.resize(27);
   for(unsigned n=0;n<27;++n) {auto p=hex.GetXcoarse(n); c.x[n]={(i+.5*(p[0]+1))/4,(j+.5*(p[1]+1))/4,(k+.5*(p[2]+1))/4};
    if(geometry=="C") {auto& x=c.x[n]; x[0]+=.4*x[0]*(1-x[0])*x[1]*(1-x[1]);}
   } cells.push_back(c);
  }
 } return cells;
}
std::vector<Cell> refine(const std::vector<Cell>& cells) {
 std::vector<Cell> out;
 for(const auto& c:cells) for(unsigned k=0;k<2;++k)for(unsigned j=0;j<2;++j)for(unsigned i=0;i<2;++i){Cell v;v.id=c.id;v.child=8*c.child+4*k+2*j+i;v.x.resize(27);for(unsigned n=0;n<27;++n){auto p=hex.GetXcoarse(n);v.x[n]=map(c,{.5*p[0]+i-.5,.5*p[1]+j-.5,.5*p[2]+k-.5});}out.push_back(v);}return out;
}
std::vector<double> nodal(const Cell& c){std::vector<double> g;for(auto x:c.x)g.push_back(sphere(x));return g;}
double legacy(const std::vector<double>& g){for(double v:g)if(g[0]*v<=0)return .5;return g[0]>0?1:0;}
void survey(std::vector<Cell> cells,std::string out) {
 std::ofstream summary(out+"/cut_counts.csv"), details(out+"/new_cut_cells.csv");
 summary<<"level,elements,nodal,bernstein,dense,missed_by_bernstein,lost_nodal\n";
 details<<"level,neu_parent_id,child_path,x,y,z,sample_x,sample_y,sample_z,sample_G,nodal_min,nodal_max\n";
 for(unsigned level=0;level<3;++level){unsigned old=0,now=0,dense=0,miss=0,lost=0;
  for(auto& c:cells){auto g=nodal(c);double a=legacy(g),b=(FELevelSetMoments::NodalLevelSetPhase(HEX,g)==0 ? .5 : (FELevelSetMoments::NodalLevelSetPhase(HEX,g)>0 ? 1. : 0.));old+=a==.5;now+=b==.5;lost+=a==.5&&b!=.5;
   // Independent direct physical-sphere evaluation: valid here only after Q1 audit.
   bool pos=false,neg=false;Point witness={};double best=-1e9;
   for(unsigned i=0;i<=12;++i)for(unsigned j=0;j<=12;++j)for(unsigned k=0;k<=12;++k){Point p={-1+i/6.,-1+j/6.,-1+k/6.};auto x=map(c,p);double v=sphere(x);pos|=v>0;neg|=v<0;if(v>best){best=v;witness=x;}}
   bool hit=pos&&neg;dense+=hit;miss+=hit&&b!=.5;
   if(b==.5&&a!=.5){auto x=map(c,{0,0,0});details<<std::setprecision(17)<<level<<','<<c.id<<','<<c.child;for(double q:x)details<<','<<q;for(double q:witness)details<<','<<q;details<<','<<best<<','<<*std::min_element(g.begin(),g.end())<<','<<*std::max_element(g.begin(),g.end())<<'\n';}
  }
  summary<<level<<','<<cells.size()<<','<<old<<','<<now<<','<<dense<<','<<miss<<','<<lost<<std::endl;
  std::cout<<"survey "<<level<<' '<<old<<' '<<now<<' '<<dense<<std::endl;
  if(miss||lost)throw std::runtime_error("Classifier missed independent crossing");
  if(level<2)cells=refine(cells);
 }
}
void audit(const std::vector<Cell>& cells) {
 double q1=0,variation=0,levelError=0,minDet=1e9;
 for(auto& c:cells){Eigen::Matrix3d first;bool init=false;
  for(unsigned n=0;n<27;++n){auto r=hex.GetXcoarse(n);Point p={r[0],r[1],r[2]},linear={};for(unsigned v=0;v<8;++v){auto t=hex.GetXcoarse(v);double w=(1+t[0]*p[0])*(1+t[1]*p[1])*(1+t[2]*p[2])/8;for(unsigned d=0;d<3;++d)linear[d]+=w*c.x[v][d];}for(unsigned d=0;d<3;++d)q1=std::max(q1,std::abs(c.x[n][d]-linear[d]));}
  auto g=nodal(c);
  for(unsigned i=0;i<3;++i)for(unsigned j=0;j<3;++j)for(unsigned k=0;k<3;++k){Point p={-.73+.73*i,-.73+.73*j,-.73+.73*k};Eigen::Matrix3d J=Eigen::Matrix3d::Zero();double gh=0;for(unsigned n=0;n<27;++n){auto d=grad(n,p);gh+=phi(n,p)*g[n];for(unsigned a=0;a<3;++a)for(unsigned b=0;b<3;++b)J(a,b)+=c.x[n][a]*d[b];}if(!init){first=J;init=true;}variation=std::max(variation,(J-first).norm());minDet=std::min(minDet,J.determinant());levelError=std::max(levelError,std::abs(gh-sphere(map(c,p))));}
 }
 std::cout<<std::setprecision(17)<<"audit q1_node_defect="<<q1<<" jacobian_variation="<<variation<<" min_det_sample="<<minDet<<" sphere_pullback_defect="<<levelError<<std::endl;
}
std::unique_ptr<CutFemWeight_LS<long double,boost::multiprecision::cpp_bin_float_quad>> cachedRule;
unsigned cachedOrder=0;
int runConfig(int argc,char**argv) {try{
 if(argc<4 || (argc<8 && std::string(argv[1])=="run"))throw std::runtime_error("usage: static_balance survey|audit|run geometry mesh [m tol legacy output_prefix]");
 auto cells=mesh(argv[2],argv[3]);
 if(std::string(argv[1])=="survey"){survey(cells,argc>4?argv[4]:".");return 0;}
 if(std::string(argv[1])=="audit"){audit(cells);return 0;}
 unsigned m=std::stoul(argv[4]);double tol=std::stod(argv[5]);bool old=std::stoi(argv[6]);std::string prefix=argv[7];
 if(m<1||m>18||tol<=0)throw std::runtime_error("Unsupported order/tolerance");
 auto start=std::chrono::steady_clock::now();
 if(!cachedRule || cachedOrder!=m){cachedRule.reset(new CutFemWeight_LS<long double,boost::multiprecision::cpp_bin_float_quad>(HEX,m,"legendre"));cachedOrder=m;}
 auto& rule=*cachedRule;rule.ClearCache();rule.SetDirectionValidation(0);
 FELevelSetMoments::OuterQuadratureRules rules;
 unsigned gk=15;
 if(const char* requested=std::getenv("CUTFEM_GK"))gk=std::stoul(requested);
 rules.regular2D=rules.tangent2D=rules.regular3D=rules.tangent3D=gk;
 rule.SetOuterQuadratureRules(rules);
 FELevelSetMoments::OuterIntervalQuadrature::Statistics(2)={};
 FELevelSetMoments::OuterIntervalQuadrature::Statistics(3)={};
 std::map<std::array<long long,3>,unsigned> ids;std::vector<Point> nodes;std::vector<std::array<unsigned,27>> dofs;
 for(auto& c:cells){std::array<unsigned,27> v;for(unsigned n=0;n<27;++n){std::array<long long,3> key;for(unsigned d=0;d<3;++d)key[d]=std::llround(c.x[n][d]*1e10);auto it=ids.emplace(key,ids.size());if(it.second)nodes.push_back(c.x[n]);v[n]=it.first->second;}dofs.push_back(v);}
 std::vector<double> S(3*ids.size()),P(S.size());double volume=0,area=0,momentTime=0,maxWeight=0;unsigned cuts=0;
 std::ofstream local(prefix+"_elements.csv"),forces(prefix+"_forces.csv");local<<"element,x,y,z,surface_l2,pressure_l2,imbalance_l2\n";
 CutFemWeight_LS<double,long double> fullRule(HEX,6,"legendre");
 unsigned nf=fullRule.GetGaussQuadraturePointNumber();
 std::vector<std::vector<double>> gf(nf,std::vector<double>(81)),pf(nf,std::vector<double>(8));
 auto pressureShape=[&](unsigned v,const Point& p){const double* t=hex.GetXcoarse(v);return (1+t[0]*p[0])*(1+t[1]*p[1])*(1+t[2]*p[2])/8.;};
 for(unsigned q=0;q<nf;++q){Point p;for(unsigned a=0;a<3;++a)p[a]=fullRule.GetGaussCoordinatePointer(a)[q];for(unsigned v=0;v<8;++v)pf[q][v]=pressureShape(v,p);for(unsigned n=0;n<27;++n){auto z=grad(n,p);for(unsigned a=0;a<3;++a)gf[q][3*n+a]=z[a];}}
 std::ofstream bt(prefix+"_bt.csv");bt<<"row,col,value\n"<<std::setprecision(17);
 std::ofstream meta(prefix+"_nodes.csv");meta<<"node,x,y,z\n"<<std::setprecision(17);for(unsigned n=0;n<nodes.size();++n)meta<<n<<','<<nodes[n][0]<<','<<nodes[n][1]<<','<<nodes[n][2]<<'\n';
 auto outputBlock=[&](unsigned e,const std::vector<double>& in,const std::vector<double>& all){for(unsigned n=0;n<27;++n)for(unsigned a=0;a<3;++a){if(nodes[dofs[e][n]][a]<1e-9||nodes[dofs[e][n]][a]>1-1e-9)continue;for(unsigned v=0;v<8;++v){unsigned j=(3*n+a)*8+v;unsigned row=3*dofs[e][n]+a,col=dofs[e][v];if(in[j]!=0)bt<<row<<','<<col<<','<<in[j]<<'\n';double out=all[j]-in[j];if(out!=0)bt<<row<<','<<nodes.size()+col<<','<<out<<'\n';}}};
 // Cache reference shapes, not geometry or moments.
 unsigned ng=rule.GetGaussQuadraturePointNumber();std::vector<std::vector<double>> gradients(ng,std::vector<double>(81));
 for(unsigned q=0;q<ng;++q){Point p;for(unsigned d=0;d<3;++d)p[d]=rule.GetGaussCoordinatePointer(d)[q];for(unsigned n=0;n<27;++n){auto d=grad(n,p);for(unsigned a=0;a<3;++a)gradients[q][3*n+a]=d[a];}}
 std::vector<std::vector<double>> ps(ng,std::vector<double>(8));for(unsigned q=0;q<ng;++q){Point p;for(unsigned a=0;a<3;++a)p[a]=rule.GetGaussCoordinatePointer(a)[q];for(unsigned v=0;v<8;++v)ps[q][v]=pressureShape(v,p);}
 for(unsigned e=0;e<cells.size();++e){auto& c=cells[e];auto g=nodal(c);double color=old?legacy(g):(FELevelSetMoments::NodalLevelSetPhase(HEX,g)==0 ? .5 : (FELevelSetMoments::NodalLevelSetPhase(HEX,g)>0 ? 1. : 0.));std::vector<double> full(648,0.),inner(648,0.);
  for(unsigned q=0;q<nf;++q){auto& d=gf[q];Eigen::Matrix3d J=Eigen::Matrix3d::Zero();for(unsigned n=0;n<27;++n)for(unsigned a=0;a<3;++a)for(unsigned b=0;b<3;++b)J(a,b)+=c.x[n][a]*d[3*n+b];double det=J.determinant();if(det<=0)throw std::runtime_error("invalid geometry");Eigen::Matrix3d invT=J.inverse().transpose();double qw=fullRule.GetGaussWeightPointer()[q];for(unsigned n=0;n<27;++n)for(unsigned a=0;a<3;++a){double phys=0;for(unsigned b=0;b<3;++b)phys+=invT(a,b)*d[3*n+b];double t=qw*det*phys;for(unsigned v=0;v<8;++v)full[(3*n+a)*8+v]+=t*pf[q][v];}}
  if(color==0){outputBlock(e,inner,full);continue;}
  std::vector<std::vector<double>> x(3,std::vector<double>(27));for(unsigned d=0;d<3;++d)for(unsigned n=0;n<27;++n)x[d][n]=c.x[n][d];
  std::vector<long double> w(ng,color),s(ng,0); std::vector<std::vector<long double>> normal; bool three=std::getenv("CUTFEM_THREE") != nullptr; std::vector<double> ls(81),lp(81);
  if(color==.5){++cuts;std::cerr<<"cut element "<<c.id<<std::endl;auto warningBefore=FELevelSetMoments::QuadratureWarningCount();auto begin=std::chrono::steady_clock::now();try{if(three) rule.GetNormalSurfaceWeightsFromNodalGeometry(g,x,w,normal,tol,tol); else rule.GetTangentialSurfaceWeightsFromNodalGeometry(g,x,w,s,tol,tol);}catch(const std::exception& er){std::cerr<<"element="<<c.id<<" "<<er.what()<<std::endl;throw;}if(FELevelSetMoments::QuadratureWarningCount()!=warningBefore) throw std::runtime_error("Moment quadrature failed its tolerance; see warning");momentTime+=std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count();}
  for(unsigned q=0;q<ng;++q){if(!std::isfinite(w[q])||!std::isfinite(s[q]))throw std::runtime_error("nonfinite weight");maxWeight=std::max(maxWeight,static_cast<double>(std::max(std::abs(w[q]),std::abs(s[q]))));
   auto& d=gradients[q];Eigen::Matrix3d J=Eigen::Matrix3d::Zero();std::vector<double> dg(3);for(unsigned n=0;n<27;++n)for(unsigned a=0;a<3;++a){dg[a]+=g[n]*d[3*n+a];for(unsigned b=0;b<3;++b)J(a,b)+=x[a][n]*d[3*n+b];}
   double det=J.determinant();if(det<=0)throw std::runtime_error("invalid geometry");Eigen::Matrix3d invT=J.inverse().transpose();std::vector<std::vector<double>> jac(3,std::vector<double>(3));for(unsigned a=0;a<3;++a)for(unsigned b=0;b<3;++b)jac[a][b]=J(a,b);
   std::vector<double> force,tangent;double metric=FELevelSetMoments::TangentialSurfaceFactors(jac,dg,d,force,tangent);double qw=rule.GetGaussWeightPointer()[q];volume+=qw*w[q]*det;area+=qw*s[q]*metric;
   for(unsigned n=0;n<27;++n)for(unsigned a=0;a<3;++a){double physical=0;for(unsigned b=0;b<3;++b)physical+=invT(a,b)*d[3*n+b];lp[3*n+a]+=10*qw*w[q]*det*physical;if(color==.5){double z=qw*w[q]*det*physical;for(unsigned v=0;v<8;++v)inner[(3*n+a)*8+v]+=z*ps[q][v];}ls[3*n+a]-=three && color==.5 ? qw*(normal[(a+1)%3][q]*tangent[3*n+(a+2)%3]-normal[(a+2)%3][q]*tangent[3*n+(a+1)%3]) : qw*s[q]*force[3*n+a];}
  }
  if(color==1)inner=full;outputBlock(e,inner,full);
  double sn=0,pn=0,bn=0;for(unsigned n=0;n<27;++n)for(unsigned d=0;d<3;++d){unsigned i=3*n+d,j=3*dofs[e][n]+d;S[j]+=ls[i];P[j]+=lp[i];sn+=ls[i]*ls[i];pn+=lp[i]*lp[i];bn+=(ls[i]+lp[i])*(ls[i]+lp[i]);}
  if(color==.5){auto xc=map(c,{0,0,0});local<<std::setprecision(17)<<c.id;for(double a:xc)local<<','<<a;local<<','<<sqrt(sn)<<','<<sqrt(pn)<<','<<sqrt(bn)<<'\n';}
 }
 double sn=0,pn=0,bn=0;forces<<"node,component,surface,pressure,imbalance\n";
 for(unsigned n=0;n<nodes.size();++n)for(unsigned d=0;d<3;++d){if(nodes[n][d]<1e-9||nodes[n][d]>1-1e-9)continue;unsigned j=3*n+d;sn+=S[j]*S[j];pn+=P[j]*P[j];bn+=(S[j]+P[j])*(S[j]+P[j]);forces<<std::setprecision(17)<<n<<','<<d<<','<<S[j]<<','<<P[j]<<','<<S[j]+P[j]<<'\n';}
 if(std::getenv("CUTFEM_THREE")) area=std::numeric_limits<double>::quiet_NaN();
 double exactV=4*M_PI*.008/3,exactA=4*M_PI*.04;
 std::ofstream csv(prefix+".csv");csv<<"mesh,refinement,m,abs_tol,rel_tol,classifier,cut_count,surface_norm,pressure_norm,imbalance_norm,balance_error,other_norm,total_norm,volume,volume_rel_error,area,area_rel_error,runtime_s,moment_s,max_abs_weight,spurious_velocity,status,gk,outer_samples_2d,outer_samples_3d\n";
 csv<<std::setprecision(17)<<argv[2]<<",0,"<<m<<','<<tol<<','<<tol<<','<<(old?"nodal":"bernstein")<<','<<cuts<<','<<sqrt(sn)<<','<<sqrt(pn)<<','<<sqrt(bn)<<','<<sqrt(bn/std::max(sn,pn))<<",0,"<<sqrt(bn)<<','<<volume<<','<<std::abs(volume-exactV)/exactV<<','<<area<<','<<std::abs(area-exactA)/exactA<<','<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<','<<momentTime<<','<<maxWeight<<",,native_assembled,"<<gk<<','<<FELevelSetMoments::OuterIntervalQuadrature::Statistics(2).samples<<','<<FELevelSetMoments::OuterIntervalQuadrature::Statistics(3).samples<<'\n';
 if(std::string(argv[1])=="test" && (std::string(argv[2])!="B" || cuts!=54 || sqrt(bn/std::max(sn,pn))>.16 || std::abs(volume-exactV)/exactV>1e-8 || std::abs(area-exactA)/exactA>2e-4))throw std::runtime_error("Unstructured sphere force-balance regression");
 std::cout<<"balance="<<sqrt(bn/std::max(sn,pn))<<" cuts="<<cuts<<std::endl;
 return 0;
 }catch(const std::exception& e){std::cerr<<e.what()<<std::endl;return 2;}}

int main(int argc,char**argv){
 if(argc>1 && std::string(argv[1])=="--batch") {
  if(argc!=4)return 2;
  std::ifstream input(argv[3]);std::string g,m,t,old,p;
  while(input>>g>>m>>t>>old>>p){
   std::vector<std::string> args={argv[0],"run",g,argv[2],m,t,old,p};std::vector<char*> v;for(auto& a:args)v.push_back(&a[0]);
   std::ofstream log(p+".log");auto* co=std::cout.rdbuf(log.rdbuf());auto* ce=std::cerr.rdbuf(log.rdbuf());
   int result=runConfig(v.size(),v.data());std::cout.rdbuf(co);std::cerr.rdbuf(ce);
   std::ofstream(p+".exit")<<result;std::cout<<p<<" exit="<<result<<std::endl;
  }return 0;
 }
 return runConfig(argc,argv);
}
