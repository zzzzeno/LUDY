#include <iostream>
#include <tgmath.h>
#include <algorithm>  
#include <complex.h>
//#include <quadmath.h>
#include <vector>
#include <stdlib.h>
#include "dual"
#include "mp++/mp++.hpp"



using namespace std;
using namespace duals::literals;

int DEBUG=0;
typedef mppp::complex128 dcomp;
typedef mppp::real128 dreal;

struct jet_list{
    vector<vector<dreal>> jet_momenta;
    vector<vector<int>> jet_labels;
};


dreal NLO_constant(dreal m){return 0.389379*pow(10,9.)*pow(4*2*acos(0.0)/132.507,2.)*(4*2*acos(0.0)*0.118/9)/(pow(2*2*acos(0.0),5.)*4*8*pow(m,6.));};

dreal norm(vector<dreal> v){return pow(pow(v[0],2.)+pow(v[1],2.)+pow(v[2],2.),0.5);};
dcomp norm(vector<dcomp> v){return pow(pow(v[0],2.)+pow(v[1],2.)+pow(v[2],2.),0.5);};
dreal spatial_norm(vector<dreal> v){return pow(pow(v[1],2.)+pow(v[2],2.)+pow(v[3],2.),0.5);};

dreal KroneckerDelta(int i, int j){if(i==j){return 1;}; return 0;};

dcomp Power(dcomp x, dreal c){return pow(x,c);}
dcomp Sqrt(dcomp x){return pow(x,0.5);}

dreal Power(dreal x, dreal c){return pow(x,c);}
dreal Sqrt(dreal x){return pow(x,0.5);}
dreal Abs(dreal x){return abs(x);}
dreal SP(vector<dreal> v1, vector<dreal> v2){return v1[0]*v2[0]-v1[1]*v2[1]-v1[2]*v2[2]-v1[3]*v2[3];};
dcomp SP(vector<dcomp> v1, vector<dcomp> v2){return v1[0]*v2[0]-v1[1]*v2[1]-v1[2]*v2[2]-v1[3]*v2[3];};
dcomp SP(vector<dcomp> v1, vector<dreal> v2){return v1[0]*v2[0]-v1[1]*v2[1]-v1[2]*v2[2]-v1[3]*v2[3];};
dcomp SP(vector<dreal> v1, vector<dcomp> v2){return v1[0]*v2[0]-v1[1]*v2[1]-v1[2]*v2[2]-v1[3]*v2[3];};

vector<dreal> vector_minus(vector<dreal> a, vector<dreal> b){
    vector<dreal> res=a;
    for(int i=0; i<a.size(); i++){res[i]=a[i]-b[i];};
    return res;
};

vector<dcomp> vector_minus(vector<dcomp> a, vector<dreal> b){
    vector<dcomp> res=a;
    for(int i=0; i<a.size(); i++){res[i]=a[i]-b[i];};
    return res;
};

vector<dreal> vector_plus(vector<dreal> a, vector<dreal> b){
    vector<dreal> res=a;
    for(int i=0; i<a.size(); i++){res[i]=a[i]+b[i];};
    return res;
};

vector<dcomp> vector_plus(vector<dcomp> a, vector<dreal> b){
    vector<dcomp> res=a;
    for(int i=0; i<a.size(); i++){res[i]=a[i]+b[i];};
    return res;
};

vector<dreal> vector_swap(vector<dreal> a){
    vector<dreal> res=a;
    for(int i=0; i<a.size(); i++){res[i]=-a[i];};
    return res;
};

vector<dreal> vector_times(vector<dreal> a, dreal c){
    vector<dreal> res=a;
    for(int i=0; i<a.size(); i++){res[i]=a[i]*c;};
    return res;    
};



struct observable{
    dreal eval;
    dreal jac;
    vector<dreal> j1;
    int spin1;
    vector<dreal> j2;
    int spin2;
    vector<dreal> pg;
};

struct observable_c{
    dcomp eval;
    dcomp jac;
    vector<dreal> j1;
    int spin1;
    vector<dreal> j2;
    int spin2;
    vector<dreal> pg;
};



class causal_flow {

    protected:
        dreal m;
        dreal sigma;

    public:
        causal_flow(){};

        causal_flow(dreal nm, dreal nsigma){m=nm; sigma=nsigma;};

        void set_parameters(dreal nm, dreal nsigma){m=nm; sigma=nsigma;};

        dreal t_val(vector<vector<dreal>> vi, vector<vector<dreal>> vf, vector<dreal> q){
            dreal ni=0;
            dreal nf=0;
            for(int i=0; i<vi.size(); i++){               
                ni=ni+norm(vi[i]);
            };

            for(int f=0; f<vf.size(); f++){
                nf=nf+norm(vf[f]);
            };

            if(pow(ni-nf,2.)-pow(norm(q),2.)<0){
                return 999;
            };

            return sqrt(pow(m,2.)/(pow(ni-nf,2.)-pow(norm(q),2.)));
        };

        dreal t_val_jj(vector<vector<dreal>> vi, dreal q0){
            dreal ni=0;
            
            for(int i=0; i<vi.size(); i++){               
                ni=ni+norm(vi[i]);
            };

            return abs(q0/ni);
        };

        vector<dreal> perform_flow(vector<vector<dreal>> vi, vector<vector<dreal>> vf, vector<dreal> q, vector<dreal> v){
            dreal t=t_val(vi, vf, q);
            vector<dreal> res;
            for (int i=0; i<v.size(); i++){
                res.push_back(t*v[i]);
            };
            return res;
        };

        vector<dreal> perform_flow_jj(vector<vector<dreal>> vi, dreal q0, vector<dreal> v){
            dreal t=t_val_jj(vi, q0);
            vector<dreal> res;
            for (int i=0; i<v.size(); i++){
                res.push_back(t*v[i]);
            };
            return res;
        };

        dreal jacques(vector<vector<dreal>> vi, vector<vector<dreal>> vf, vector<dreal> q, int n_loops){

            dreal t=t_val(vi,vf,q);

            dreal ni=0;
            dreal nf=0;
            for(int i=0; i<vi.size(); i++){
                ni+=norm(vi[i]);
            };
            for(int f=0; f<vf.size(); f++){
                nf+=norm(vf[f]);
            };

            return pow(t,3*n_loops)/abs(ni-nf-t*pow(norm(q),2.)/sqrt(pow(t,2.)*pow(norm(q),2.)+pow(m,2)));

        };

        dreal jacques_jj(vector<vector<dreal>> vi, dreal q0, int n_loops){

            dreal t=t_val_jj(vi,q0);

            dreal ni=0;

            for(int i=0; i<vi.size(); i++){
                ni+=norm(vi[i]);
            };

            return pow(t,3*n_loops)*q0/ni;

        };

        dreal h(dreal t){
            return 2/(sqrt(2*M_PI)*sigma)*exp(-pow(t,2.)/(2*pow(sigma,2.)));
        };

        /*dreal h(dreal t){
            dreal n=0.11993777196806144736803650163679;
            return 1/n * exp(-pow(t,2.)-1/pow(t,2.));
        };*/

};


class boost_flow {

    private:

    public:
        boost_flow(){};

        dreal t_val(vector<dreal> j1pj2, int comp){
            dreal p0=j1pj2[0];
            dreal pi=j1pj2[comp];

            return pi/p0;
        };

        vector<dreal> perform_flow(vector<dreal> j1pj2, int comp, dreal t){
            dreal x=t;

            dreal gamma=1/sqrt(1-pow(x,2.));

            vector<vector<dreal>> boost={{gamma,-KroneckerDelta(comp,1)*x*gamma,-KroneckerDelta(comp,2)*x*gamma,-KroneckerDelta(comp,3)*x*gamma},
            {-KroneckerDelta(comp,1)*x*gamma,1+(gamma-1)*KroneckerDelta(comp,1),0,0},
            {-KroneckerDelta(comp,2)*x*gamma,0,1+(gamma-1)*KroneckerDelta(comp,2),0},
            {-KroneckerDelta(comp,3)*x*gamma,0,0,1+(gamma-1)*KroneckerDelta(comp,3)}};

            vector<dreal> res={0,0,0,0};

            for(int i=0; i<4; i++){
                dreal compi=0;
                for(int j=0; j<4; j++){
                    compi+=boost[i][j]*j1pj2[j];
                };    
                res[i]=compi;
            };

            return res;
        };

        vector<dreal> perform_jacques(vector<dreal> v, dreal x, int comp){

            dreal gamma=1/sqrt(1-pow(x,2.));
            dreal gammad=x/pow(sqrt(1-pow(x,2)),3);

            vector<vector<dreal>> boost_j={{gammad,-KroneckerDelta(comp,1)*(x*gammad+gamma),-KroneckerDelta(comp,2)*(x*gammad+gamma),-KroneckerDelta(comp,3)*(x*gammad+gamma)},
            {-KroneckerDelta(comp,1)*(x*gammad+gamma),gammad*KroneckerDelta(comp,1),0,0},
            {-KroneckerDelta(comp,2)*(x*gammad+gamma),0,gammad*KroneckerDelta(comp,2),0},
            {-KroneckerDelta(comp,3)*(x*gammad+gamma),0,0,gammad*KroneckerDelta(comp,3)}};

            vector<dreal> res={0,0,0,0};

            for(int i=0; i<4; i++){
                dreal compi=0;
                for(int j=0; j<4; j++){
                    compi+=boost_j[i][j]*v[j];
                };    
                res[i]=compi;
            };

            return res;

        };

        dreal h(dreal t){
            dreal nconst=0.34029423827512622624524851744354;
            return (1/nconst)*exp(-1/(pow(pow(t,2.) - 1,2.)));
        };

        /*dreal h(dreal t){
            dreal nconst=0.0004980185228640048;
            return (1/nconst)*exp(-1/(pow(t*(pow(t,2.) - 1),2.)));
        };*/

};


class rotation_flow {

    private:

    public:
        rotation_flow(){};

        dreal t_val(vector<dreal> j1mj2, int comp){
            dreal x=j1mj2[comp];
            dreal y=j1mj2[comp+1];

            return atan2(x,y);
        };

        vector<dreal> perform_flow(vector<dreal> j1mj2, int comp, dreal t){
            dreal th=t;

            vector<vector<dreal>> rotation={{1,0,0,0},
            {0,KroneckerDelta(comp,2)+KroneckerDelta(comp,1)*cos(th),-KroneckerDelta(comp,1)*sin(th),0},
            {0,KroneckerDelta(comp,1)*sin(th),cos(th),-KroneckerDelta(comp,2)*sin(th)},
            {0,0,KroneckerDelta(comp,2)*sin(th),KroneckerDelta(comp,1)+KroneckerDelta(comp,2)*cos(th)}};
           

            vector<dreal> res={0,0,0,0};

            for(int i=0; i<4; i++){
                dreal compi=0;
                for(int j=0; j<4; j++){
                    compi+=rotation[i][j]*j1mj2[j];
                };    
                res[i]=compi;
            };

            return res;
        };

        vector<dreal> perform_jacques(vector<dreal> v, dreal th, int comp){

            vector<vector<dreal>> rotation_j={{0,0,0,0},
            {0,-KroneckerDelta(comp,1)*sin(th),-KroneckerDelta(comp,1)*cos(th),0},
            {0,KroneckerDelta(comp,1)*cos(th),-sin(th),-KroneckerDelta(comp,2)*cos(th)},
            {0,0,KroneckerDelta(comp,2)*cos(th),-KroneckerDelta(comp,2)*sin(th)}};
           
            vector<dreal> res={0,0,0,0};

            for(int i=0; i<4; i++){
                dreal compi=0;
                for(int j=0; j<4; j++){
                    compi+=rotation_j[i][j]*(v[j]);
                };    
                res[i]=compi;
            };

            return res;

        };

        dreal h(dreal th){
            dreal nconst=1.5707963267948987784450309845852;
            return (1/(2*nconst))*pow(sin(2*th),2);
            //return (1/(2*nconst))*pow(sin(2*th),4);
            //return exp(-1/pow(sin(2*th),2));
            //dreal nconst=4;
            //return (1/(2*nconst))*abs(sin(2*th));
        };

};






class observables{

    private:
        dreal eCM;
        dreal m;
        
        dreal res_c;
        dreal res_s;

        boost_flow b_flow;
        rotation_flow r_flow;

    public:
        observables(){};
        observables(dreal nres_c, dreal nres_s, dreal neCM, dreal nM) : b_flow(), r_flow() {
            res_c=nres_c;
            res_s=nres_s;
            m=nM;
            eCM=neCM;
        }

        void set_res(dreal nres_c, dreal nres_s){
            res_c=nres_c;
            res_s=nres_s;
        }

        void set_kinematics(dreal nm, dreal neCM){eCM=neCM; m=nm;};
        
        dreal SP(vector<dreal> w, vector<dreal> v){
            return v[0]*w[0]-v[1]*w[1]-v[2]*w[2]-v[3]*w[3];
        }

        dreal distance(vector<dreal> vec1, vector<dreal> vec2){
            dreal d=SP(vec1,vec2)/pow(eCM,2);///(vec1[0]*vec2[0]);
            return d;
        };


        jet_list two_jet_clustering(vector<vector<dreal>> constituents, int n_constituents){
            
            jet_list jets;

            vector<vector<int>> labels;

            for(int i=0; i<n_constituents; i++){
                labels.push_back({i});
            };

            int inter_n=n_constituents;

            vector<vector<dreal>> inter_constituents=constituents;

            while(inter_n>0){

                if(DEBUG==1){cout<<"inter_n  "<<inter_n<<endl;};

                int label1=0;
                int label2=0;
                dreal current_d=-1;

                if(inter_n==2){
                    jets.jet_momenta.push_back(inter_constituents[0]);
                    jets.jet_labels.push_back(labels[0]);
                    jets.jet_momenta.push_back(inter_constituents[1]);
                    jets.jet_labels.push_back(labels[1]);

                    return jets;
                }else{

                    for(int j1=0; j1<inter_n; j1++){
                        for(int j2=j1+1; j2<inter_n; j2++){
                            
                                dreal d=distance(inter_constituents[j1],inter_constituents[j2]);
                                if(DEBUG==1){cout<<j1<<"  "<<j2<<"  -  distance: "<<d<<endl;};
                                if(d<current_d || current_d<0){
                                    current_d=d;
                                    label1=j1;
                                    label2=j2;
                                };

                            };
                        };

                    if(DEBUG==1){cout<<"min_distance: "<<current_d<<"  R: "<<res_c<<endl;
                                 cout<<"min_distance_labels: "<<label1<<" "<<label2<<endl;};  

                    if(current_d<res_c){
                        if(DEBUG==1){cout<<"entered"<<endl;};  
                        vector<dreal> constituent_sum={inter_constituents[label1][0]+inter_constituents[label2][0],
                            inter_constituents[label1][1]+inter_constituents[label2][1],
                            inter_constituents[label1][2]+inter_constituents[label2][2],
                            inter_constituents[label1][3]+inter_constituents[label2][3]};
                        if(label1<label2){
                            inter_constituents.erase(inter_constituents.begin()+label1);
                            inter_constituents.erase(inter_constituents.begin()+(label2-1));
                        }else{
                            inter_constituents.erase(inter_constituents.begin()+label2);
                            inter_constituents.erase(inter_constituents.begin()+(label1-1));
                        };

                        inter_constituents.push_back(constituent_sum);

                        vector<int> union_labels;
                        union_labels.insert( union_labels.end(), labels[label1].begin(), labels[label1].end() );
                        union_labels.insert( union_labels.end(), labels[label2].begin(), labels[label2].end() );

                        if(label1<label2){
                            labels.erase(labels.begin()+label1);
                            labels.erase(labels.begin()+(label2-1));
                        }else{
                            labels.erase(labels.begin()+label2);
                            labels.erase(labels.begin()+(label1-1));
                        };

                        labels.push_back(union_labels);
                        inter_n--;

                    }else{
                        for(int j=0; j<inter_constituents.size(); j++){
                            jets.jet_momenta.push_back(inter_constituents[j]);
                            jets.jet_labels.push_back(labels[j]);
                        };

                        return jets;
                    }

                };

            };
            
            return jets;

        };


        vector<dreal> b2b_sampling(vector<vector<dreal>> constituents, int n_constituents, int* spins, dreal a){


            jet_list my_jets=two_jet_clustering(constituents, n_constituents);

            if(DEBUG==1){
                for(int n=0; n<my_jets.jet_momenta.size(); n++){
                    cout<<"jet number:  "<<n<<endl;
                    for(int i=0; i<4; i++){
                        cout<<my_jets.jet_momenta[n][i]<<" ";
                    };
                    cout<<endl;
                };
            };

            int size=0;

            vector<int> jet_n;

            for(int i=0; i<my_jets.jet_labels.size(); i++){
                if(my_jets.jet_momenta[i][0]>res_s){
                    size++;
                    jet_n.push_back(i);
                };
            };

            if(DEBUG==1){
                cout<<"number of jets"<<endl;
                cout<<my_jets.jet_labels.size()<<endl;
                cout<<"number of hard jets"<<endl;
                cout<<size<<endl;
            };


            if(size==2){

                dreal tot_spin1=0;
                dreal tot_spin2=0;

                for(int j1=0; j1<my_jets.jet_labels[jet_n[0]].size(); j1++){
                    tot_spin1=tot_spin1+spins[my_jets.jet_labels[jet_n[0]][j1]];
                };

                for(int j1=0; j1<my_jets.jet_labels[jet_n[1]].size(); j1++){
                    tot_spin2+=tot_spin2+spins[my_jets.jet_labels[jet_n[1]][j1]];
                };

                vector<dreal> jet1=my_jets.jet_momenta[jet_n[0]];
                vector<dreal> jet2=my_jets.jet_momenta[jet_n[1]];

                vector<dreal> j1pj2={0,0,0,0};
                vector<dreal> j1mj2={0,0,0,0};

                for(int i=0; i<4; i++){
                    j1pj2[i]=my_jets.jet_momenta[jet_n[0]][i]+my_jets.jet_momenta[jet_n[1]][i];
                    j1mj2[i]=my_jets.jet_momenta[jet_n[0]][i]-my_jets.jet_momenta[jet_n[1]][i];
                };

                //j1pj2[3]=0;
                //j1mj2[3]=my_jets.jet_momenta[jet_n[0]][3]-my_jets.jet_momenta[jet_n[1]][3];

                vector<dreal> interjb=j1pj2;
                vector<dreal> interjr=j1mj2;

                vector<dreal> interj_jacques;

                vector<dreal> vs={0,0,0};

                    if(DEBUG==1){
                        cout<<"j1+j2: "<<j1pj2[0]<<" "<<j1pj2[1]<<" "<<j1pj2[2]<<" "<<j1pj2[3]<<" "<<endl;
                        cout<<"j1-j2: "<<j1mj2[0]<<" "<<j1mj2[1]<<" "<<j1mj2[2]<<" "<<j1mj2[3]<<" "<<endl;
                    };

                dreal hs=1;
                vector<vector<dreal>> jacques_vecs={j1pj2,j1pj2,j1pj2,j1pj2,j1mj2,j1mj2,j1mj2};

                dreal tLz=b_flow.t_val(jacques_vecs[0],3);
                jacques_vecs={b_flow.perform_flow(jacques_vecs[0],3,tLz),b_flow.perform_flow(jacques_vecs[1],3,tLz),b_flow.perform_flow(jacques_vecs[2],3,tLz),b_flow.perform_jacques(jacques_vecs[3],tLz,3),b_flow.perform_flow(jacques_vecs[4],3,tLz),b_flow.perform_flow(jacques_vecs[5],3,tLz),b_flow.perform_flow(jacques_vecs[6],3,tLz)};


                dreal tLy=b_flow.t_val(jacques_vecs[0],2);
                jacques_vecs={b_flow.perform_flow(jacques_vecs[0],2,tLy),b_flow.perform_flow(jacques_vecs[1],2,tLy),b_flow.perform_jacques(jacques_vecs[2],tLy,2),b_flow.perform_flow(jacques_vecs[3],2,tLy),b_flow.perform_flow(jacques_vecs[4],2,tLy),b_flow.perform_flow(jacques_vecs[5],2,tLy),b_flow.perform_flow(jacques_vecs[6],2,tLy)};


                dreal tLx=b_flow.t_val(jacques_vecs[0],1);
                jacques_vecs={b_flow.perform_flow(jacques_vecs[0],1,tLx),b_flow.perform_jacques(jacques_vecs[1],tLx,1),b_flow.perform_flow(jacques_vecs[2],1,tLx),b_flow.perform_flow(jacques_vecs[3],1,tLx),b_flow.perform_flow(jacques_vecs[4],1,tLx),b_flow.perform_flow(jacques_vecs[5],1,tLx),b_flow.perform_flow(jacques_vecs[6],1,tLx)};


                dreal txy=r_flow.t_val(jacques_vecs[4],1);
                jacques_vecs={r_flow.perform_flow(jacques_vecs[0],1,txy),r_flow.perform_flow(jacques_vecs[1],1,txy),r_flow.perform_flow(jacques_vecs[2],1,txy),r_flow.perform_flow(jacques_vecs[3],1,txy),r_flow.perform_flow(jacques_vecs[4],1,txy),r_flow.perform_jacques(jacques_vecs[5],txy,1),r_flow.perform_flow(jacques_vecs[6],1,txy)};


                dreal tyz=r_flow.t_val(jacques_vecs[4],2);
                jacques_vecs={r_flow.perform_flow(jacques_vecs[0],2,tyz),r_flow.perform_flow(jacques_vecs[1],2,tyz),r_flow.perform_flow(jacques_vecs[2],2,tyz),r_flow.perform_flow(jacques_vecs[3],2,tyz),r_flow.perform_flow(jacques_vecs[4],2,tyz),r_flow.perform_flow(jacques_vecs[5],2,tyz),r_flow.perform_jacques(jacques_vecs[6],tyz,2)};


                dreal tLz_gen=-a/sqrt(1+pow(a,2.));
                jacques_vecs={b_flow.perform_flow(jacques_vecs[0],3,tLz_gen),b_flow.perform_flow(jacques_vecs[1],3,tLz_gen),b_flow.perform_flow(jacques_vecs[2],3,tLz_gen),b_flow.perform_flow(jacques_vecs[3],3,tLz_gen),b_flow.perform_flow(jacques_vecs[4],3,tLz_gen),b_flow.perform_flow(jacques_vecs[5],3,tLz_gen),b_flow.perform_flow(jacques_vecs[6],3,tLz_gen)};

                /*jacques_vecs={b_flow.perform_flow(jacques_vecs[0],3,tLz),b_flow.perform_flow(jacques_vecs[1],3,tLz),b_flow.perform_flow(jacques_vecs[2],3,tLz),b_flow.perform_flow(jacques_vecs[3],3,tLz),b_flow.perform_flow(jacques_vecs[4],3,tLz),b_flow.perform_flow(jacques_vecs[5],3,tLz)};*/

                dreal res_f=Sqrt(SP(jacques_vecs[0],jacques_vecs[0]))*(b_flow.h(tLx)*b_flow.h(tLy)*b_flow.h(tLz)*r_flow.h(txy)*r_flow.h(tyz))/abs(jacques_vecs[1][1]*jacques_vecs[2][2]*jacques_vecs[3][3]*jacques_vecs[5][1]*jacques_vecs[6][2]);

                if(DEBUG==1){
                    for(int i=0; i<7; i++){
                        cout<<"vec "<<i<< ":   ";
                        for(int j=0; j<4; j++){
                            cout<<jacques_vecs[i][j]<<" ";
                        };
                        cout<<endl;
                    };

                    cout<<"--------------------"<<endl;

                    cout<<"ts"<<endl;
                    cout<<"Lx: "<<tLx<<endl;
                    cout<<"Ly: "<<tLy<<endl;
                    cout<<"Ly: "<<tLz<<endl;
                    cout<<"Rxy: "<<txy<<endl;
                    cout<<"Ryz: "<<tyz<<endl;

                    cout<<"jacques transf"<<endl;
                    cout<<"Lx: "<<jacques_vecs[1][1]<<endl;
                    cout<<"Ly: "<<jacques_vecs[2][2]<<endl;
                    cout<<"Lz: "<<jacques_vecs[3][3]<<endl;
                    cout<<"Rxy: "<<jacques_vecs[5][1]<<endl;
                    cout<<"Ryz: "<<jacques_vecs[6][2]<<endl;

                    cout<<"hts"<<endl;
                    cout<<"Lx: "<<b_flow.h(tLx)<<endl;
                    cout<<"Ly: "<<b_flow.h(tLy)<<endl;
                    cout<<"Lz: "<<b_flow.h(tLz)<<endl;
                    cout<<"Rxy: "<<r_flow.h(txy)<<endl;
                    cout<<"Ryz: "<<r_flow.h(tyz)<<endl;

                    cout<<"jac   "<<res_f<<endl;
                };

                //dreal x1px2=jacques_vecs[0][0]/eCM;
                //dreal x1mx2=jacques_vecs[3][0]/eCM;

                dreal x1px2=(jacques_vecs[0][0])/eCM;
                dreal x1mx2=(jacques_vecs[0][3])/eCM;

                if(DEBUG==1){
                    cout<<"j1+j2 post: "<<jacques_vecs[0][0]<<" "<<jacques_vecs[0][1]<<" "<<jacques_vecs[0][2]<<" "<<jacques_vecs[0][3]<<" "<<endl;
                    cout<<"j1-j2 post: "<<jacques_vecs[4][0]<<" "<<jacques_vecs[4][1]<<" "<<jacques_vecs[4][2]<<" "<<jacques_vecs[4][3]<<" "<<endl;
                    cout<<"x1: "<<(x1mx2+x1px2)/2<<" x2: "<<(x1px2-x1mx2)/2<<endl;
                };

                vector<dreal> res_vec;


                res_vec={res_f, (x1mx2+x1px2)/2, (x1px2-x1mx2)/2, tot_spin1, tot_spin2};

                return res_vec;

            };

            vector<dreal> res_vec_n={0,0,0,0,0};

            return res_vec_n;

        };





        vector<dreal> b2b_sampling_pt(vector<vector<dreal>> constituents, int n_constituents, vector<dreal> final_parton, int* spins, dreal a){

            //dreal old_res_c=res_c;
            //dreal old_res_s=res_s;

            /*if(n_constituents==2){
                res_c=0;
                res_s=0;
            };*/

            jet_list my_jets=two_jet_clustering(constituents, n_constituents);

            /*if(n_constituents==2){
                res_c=old_res_c;
                res_s=old_res_s;
            };*/

            if(DEBUG==1){
                for(int n=0; n<my_jets.jet_momenta.size(); n++){
                    cout<<"jet number:  "<<n<<endl;
                    for(int i=0; i<4; i++){
                        cout<<my_jets.jet_momenta[n][i]<<" ";
                    };
                    cout<<endl;
                };
            };

            
            int size=0;
            vector<int> jet_n;


            //IS THIS REALLY NEEDED? SOFT PARTICLES WILL ALWAYS BE CLUSTERED ANYWAY WITH JADE
            /*
            if(my_jets.jet_labels.size()>2){
                for(int i=0; i<my_jets.jet_labels.size(); i++){
                    if(my_jets.jet_momenta[i][0]>res_s){
                    size++;
                    jet_n.push_back(i);
                    };
                };
            };*/

            for(int i=0; i<my_jets.jet_labels.size(); i++){  
                size++;
                jet_n.push_back(i);
                };


            if(DEBUG==1){
                cout<<"number of jets"<<endl;
                cout<<my_jets.jet_labels.size()<<endl;
                cout<<"number of hard jets"<<endl;
                cout<<size<<endl;
            };


            if(size==2){

                dreal tot_spin1=0;
                dreal tot_spin2=0;

                for(int j1=0; j1<my_jets.jet_labels[jet_n[0]].size(); j1++){
                    tot_spin1=tot_spin1+spins[my_jets.jet_labels[jet_n[0]][j1]];
                };

                for(int j1=0; j1<my_jets.jet_labels[jet_n[1]].size(); j1++){
                    tot_spin2+=tot_spin2+spins[my_jets.jet_labels[jet_n[1]][j1]];
                };

                vector<dreal> jet1=my_jets.jet_momenta[jet_n[0]];
                vector<dreal> jet2=my_jets.jet_momenta[jet_n[1]];

                vector<dreal> j1pj2={0,0,0,0};
                vector<dreal> j1mj2={0,0,0,0};

                for(int i=0; i<4; i++){
                    j1pj2[i]=my_jets.jet_momenta[jet_n[0]][i]+my_jets.jet_momenta[jet_n[1]][i];
                    j1mj2[i]=my_jets.jet_momenta[jet_n[0]][i]-my_jets.jet_momenta[jet_n[1]][i];
                };

                //j1pj2[3]=0;
                //j1mj2[3]=my_jets.jet_momenta[jet_n[0]][3]-my_jets.jet_momenta[jet_n[1]][3];

                vector<dreal> interjb=j1pj2;
                vector<dreal> interjr=j1mj2;

                vector<dreal> interj_jacques;

                vector<dreal> vs={0,0,0};

                    if(DEBUG==1){
                        cout<<"j1+j2: "<<j1pj2[0]<<" "<<j1pj2[1]<<" "<<j1pj2[2]<<" "<<j1pj2[3]<<" "<<endl;
                        cout<<"j1-j2: "<<j1mj2[0]<<" "<<j1mj2[1]<<" "<<j1mj2[2]<<" "<<j1mj2[3]<<" "<<endl;
                        cout<<"l pre: "<<final_parton[0]<<" "<<final_parton[1]<<" "<<final_parton[2]<<" "<<final_parton[3]<<" "<<endl;
                        cout<<"j1: "<<(j1pj2[0]+j1mj2[0])/2<<" "<<(j1pj2[1]+j1mj2[1])/2<<" "<<(j1pj2[2]+j1mj2[2])/2<<" "<<(j1pj2[3]+j1mj2[3])/2<<" "<<endl;
                        cout<<"j2: "<<(j1pj2[0]-j1mj2[0])/2<<" "<<(j1pj2[1]-j1mj2[1])/2<<" "<<(j1pj2[2]-j1mj2[2])/2<<" "<<(j1pj2[3]-j1mj2[3])/2<<" "<<endl;
                    };

                dreal hs=1;
                vector<vector<dreal>> jacques_vecs={j1pj2,j1pj2,j1pj2,j1pj2,j1mj2,j1mj2,j1mj2,final_parton};

                dreal tLz=b_flow.t_val(jacques_vecs[0],3);
                jacques_vecs={b_flow.perform_flow(jacques_vecs[0],3,tLz),b_flow.perform_flow(jacques_vecs[1],3,tLz),b_flow.perform_flow(jacques_vecs[2],3,tLz),b_flow.perform_jacques(jacques_vecs[3],tLz,3),b_flow.perform_flow(jacques_vecs[4],3,tLz),b_flow.perform_flow(jacques_vecs[5],3,tLz),b_flow.perform_flow(jacques_vecs[6],3,tLz),b_flow.perform_flow(jacques_vecs[7],3,tLz)};


                dreal tLy=b_flow.t_val(jacques_vecs[0],2);
                jacques_vecs={b_flow.perform_flow(jacques_vecs[0],2,tLy),b_flow.perform_flow(jacques_vecs[1],2,tLy),b_flow.perform_jacques(jacques_vecs[2],tLy,2),b_flow.perform_flow(jacques_vecs[3],2,tLy),b_flow.perform_flow(jacques_vecs[4],2,tLy),b_flow.perform_flow(jacques_vecs[5],2,tLy),b_flow.perform_flow(jacques_vecs[6],2,tLy),b_flow.perform_flow(jacques_vecs[7],2,tLy)};


                dreal tLx=b_flow.t_val(jacques_vecs[0],1);
                jacques_vecs={b_flow.perform_flow(jacques_vecs[0],1,tLx),b_flow.perform_jacques(jacques_vecs[1],tLx,1),b_flow.perform_flow(jacques_vecs[2],1,tLx),b_flow.perform_flow(jacques_vecs[3],1,tLx),b_flow.perform_flow(jacques_vecs[4],1,tLx),b_flow.perform_flow(jacques_vecs[5],1,tLx),b_flow.perform_flow(jacques_vecs[6],1,tLx),b_flow.perform_flow(jacques_vecs[7],1,tLx)};


                dreal txy=r_flow.t_val(jacques_vecs[4],1);
                jacques_vecs={r_flow.perform_flow(jacques_vecs[0],1,txy),r_flow.perform_flow(jacques_vecs[1],1,txy),r_flow.perform_flow(jacques_vecs[2],1,txy),r_flow.perform_flow(jacques_vecs[3],1,txy),r_flow.perform_flow(jacques_vecs[4],1,txy),r_flow.perform_jacques(jacques_vecs[5],txy,1),r_flow.perform_flow(jacques_vecs[6],1,txy),r_flow.perform_flow(jacques_vecs[7],1,txy)};


                dreal tyz=r_flow.t_val(jacques_vecs[4],2);
                jacques_vecs={r_flow.perform_flow(jacques_vecs[0],2,tyz),r_flow.perform_flow(jacques_vecs[1],2,tyz),r_flow.perform_flow(jacques_vecs[2],2,tyz),r_flow.perform_flow(jacques_vecs[3],2,tyz),r_flow.perform_flow(jacques_vecs[4],2,tyz),r_flow.perform_flow(jacques_vecs[5],2,tyz),r_flow.perform_jacques(jacques_vecs[6],tyz,2),r_flow.perform_flow(jacques_vecs[7],2,tyz)};


                dreal tLz_gen=-a/sqrt(1+pow(a,2.));
                jacques_vecs={b_flow.perform_flow(jacques_vecs[0],3,tLz_gen),b_flow.perform_flow(jacques_vecs[1],3,tLz_gen),b_flow.perform_flow(jacques_vecs[2],3,tLz_gen),b_flow.perform_flow(jacques_vecs[3],3,tLz_gen),b_flow.perform_flow(jacques_vecs[4],3,tLz_gen),b_flow.perform_flow(jacques_vecs[5],3,tLz_gen),b_flow.perform_flow(jacques_vecs[6],3,tLz_gen),b_flow.perform_flow(jacques_vecs[7],3,tLz_gen)};


                /*jacques_vecs={b_flow.perform_flow(jacques_vecs[0],3,tLz),b_flow.perform_flow(jacques_vecs[1],3,tLz),b_flow.perform_flow(jacques_vecs[2],3,tLz),b_flow.perform_flow(jacques_vecs[3],3,tLz),b_flow.perform_flow(jacques_vecs[4],3,tLz),b_flow.perform_flow(jacques_vecs[5],3,tLz)};*/

                dreal res_f=Sqrt(SP(jacques_vecs[0],jacques_vecs[0]))*(b_flow.h(tLx)*b_flow.h(tLy)*b_flow.h(tLz)*r_flow.h(txy)*r_flow.h(tyz))/abs(jacques_vecs[1][1]*jacques_vecs[2][2]*jacques_vecs[3][3]*jacques_vecs[5][1]*jacques_vecs[6][2]);

                if(DEBUG==1){
                    for(int i=0; i<7; i++){
                        cout<<"vec "<<i<< ":   ";
                        for(int j=0; j<4; j++){
                            cout<<jacques_vecs[i][j]<<" ";
                        };
                        cout<<endl;
                    };

                    cout<<"--------------------"<<endl;

                    cout<<"ts"<<endl;
                    cout<<"Lx: "<<tLx<<endl;
                    cout<<"Ly: "<<tLy<<endl;
                    cout<<"Ly: "<<tLz<<endl;
                    cout<<"Rxy: "<<txy<<endl;
                    cout<<"Ryz: "<<tyz<<endl;

                    cout<<"jacques transf"<<endl;
                    cout<<"Lx: "<<jacques_vecs[1][1]<<endl;
                    cout<<"Ly: "<<jacques_vecs[2][2]<<endl;
                    cout<<"Lz: "<<jacques_vecs[3][3]<<endl;
                    cout<<"Rxy: "<<jacques_vecs[5][1]<<endl;
                    cout<<"Ryz: "<<jacques_vecs[6][2]<<endl;

                    cout<<"hts"<<endl;
                    cout<<"Lx: "<<b_flow.h(tLx)<<endl;
                    cout<<"Ly: "<<b_flow.h(tLy)<<endl;
                    cout<<"Lz: "<<b_flow.h(tLz)<<endl;
                    cout<<"Rxy: "<<r_flow.h(txy)<<endl;
                    cout<<"Ryz: "<<r_flow.h(tyz)<<endl;

                    cout<<"jac   "<<res_f<<endl;
                };

                //dreal x1px2=jacques_vecs[0][0]/eCM;
                //dreal x1mx2=jacques_vecs[3][0]/eCM;

                dreal x1px2=2*(jacques_vecs[0][0])/eCM;
                dreal x1mx2=2*(jacques_vecs[0][3])/eCM;
                dreal pt=Sqrt(pow(jacques_vecs[7][1],2.)+pow(jacques_vecs[7][2],2.));

                if(DEBUG==1){
                    cout<<"j1+j2 post: "<<jacques_vecs[0][0]<<" "<<jacques_vecs[0][1]<<" "<<jacques_vecs[0][2]<<" "<<jacques_vecs[0][3]<<" "<<endl;
                    cout<<"j1-j2 post: "<<jacques_vecs[4][0]<<" "<<jacques_vecs[4][1]<<" "<<jacques_vecs[4][2]<<" "<<jacques_vecs[4][3]<<" "<<endl;
                    cout<<"l post: "<<jacques_vecs[7][0]<<" "<<jacques_vecs[7][1]<<" "<<jacques_vecs[7][2]<<" "<<jacques_vecs[7][3]<<" "<<endl;
                    cout<<"j1 post: "<<(jacques_vecs[0][0]+jacques_vecs[4][0])/2<<" "<<(jacques_vecs[0][1]+jacques_vecs[4][1])/2<<" "<<(jacques_vecs[0][2]+jacques_vecs[4][2])/2<<" "<<(jacques_vecs[0][3]+jacques_vecs[4][3])/2<<" "<<endl;
                    cout<<"j2 post: "<<(jacques_vecs[0][0]-jacques_vecs[4][0])/2<<" "<<(jacques_vecs[0][1]-jacques_vecs[4][1])/2<<" "<<(jacques_vecs[0][2]-jacques_vecs[4][2])/2<<" "<<(jacques_vecs[0][3]-jacques_vecs[4][3])/2<<" "<<endl;
                    cout<<"x1: "<<(x1mx2+x1px2)/2<<" x2: "<<(x1px2-x1mx2)/2<<endl;
                };

                vector<dreal> res_vec;


                res_vec={res_f, (x1mx2+x1px2)/2, (x1px2-x1mx2)/2, tot_spin1, tot_spin2, pt};

                return res_vec;

            };

            vector<dreal> res_vec_n={0,0,0,0,0,0};
            
            return res_vec_n;

        };


        observable b2b_sampling_final(vector<vector<dreal>> constituents, int n_constituents, vector<dreal> final_parton, int* spins, dreal a){

            //dreal old_res_c=res_c;
            //dreal old_res_s=res_s;

            /*if(n_constituents==2){
                res_c=0;
                res_s=0;
            };*/

            jet_list my_jets=two_jet_clustering(constituents, n_constituents);

            /*if(n_constituents==2){
                res_c=old_res_c;
                res_s=old_res_s;
            };*/

            if(DEBUG==1){
                for(int n=0; n<my_jets.jet_momenta.size(); n++){
                    cout<<"jet number:  "<<n<<endl;
                    for(int i=0; i<4; i++){
                        cout<<my_jets.jet_momenta[n][i]<<" ";
                    };
                    cout<<endl;
                };
            };

            
            int size=0;
            vector<int> jet_n;


            //IS THIS REALLY NEEDED? SOFT PARTICLES WILL ALWAYS BE CLUSTERED ANYWAY WITH JADE
            /*
            if(my_jets.jet_labels.size()>2){
                for(int i=0; i<my_jets.jet_labels.size(); i++){
                    if(my_jets.jet_momenta[i][0]>res_s){
                    size++;
                    jet_n.push_back(i);
                    };
                };
            };*/

            for(int i=0; i<my_jets.jet_labels.size(); i++){  
                size++;
                jet_n.push_back(i);
                };


            if(DEBUG==1){
                cout<<"number of jets"<<endl;
                cout<<my_jets.jet_labels.size()<<endl;
                cout<<"number of hard jets"<<endl;
                cout<<size<<endl;
            };


            if(size==2){

                dreal tot_spin1=0;
                dreal tot_spin2=0;

                for(int j1=0; j1<my_jets.jet_labels[jet_n[0]].size(); j1++){
                    tot_spin1=tot_spin1+spins[my_jets.jet_labels[jet_n[0]][j1]];
                };

                for(int j1=0; j1<my_jets.jet_labels[jet_n[1]].size(); j1++){
                    tot_spin2=tot_spin2+spins[my_jets.jet_labels[jet_n[1]][j1]];
                };


                vector<dreal> jet1=my_jets.jet_momenta[jet_n[0]];
                vector<dreal> jet2=my_jets.jet_momenta[jet_n[1]];

                vector<dreal> j1pj2={0,0,0,0};
                vector<dreal> j1mj2={0,0,0,0};

                for(int i=0; i<4; i++){
                    j1pj2[i]=my_jets.jet_momenta[jet_n[0]][i]+my_jets.jet_momenta[jet_n[1]][i];
                    j1mj2[i]=my_jets.jet_momenta[jet_n[0]][i]-my_jets.jet_momenta[jet_n[1]][i];
                };

                //j1pj2[3]=0;
                //j1mj2[3]=my_jets.jet_momenta[jet_n[0]][3]-my_jets.jet_momenta[jet_n[1]][3];

                vector<dreal> interjb=j1pj2;
                vector<dreal> interjr=j1mj2;

                vector<dreal> interj_jacques;

                vector<dreal> vs={0,0,0};

                    if(DEBUG==4){
                        cout<<"j1+j2: "<<j1pj2[0]<<" "<<j1pj2[1]<<" "<<j1pj2[2]<<" "<<j1pj2[3]<<" "<<endl;
                        cout<<"j1-j2: "<<j1mj2[0]<<" "<<j1mj2[1]<<" "<<j1mj2[2]<<" "<<j1mj2[3]<<" "<<endl;
                        cout<<"l pre: "<<final_parton[0]<<" "<<final_parton[1]<<" "<<final_parton[2]<<" "<<final_parton[3]<<" "<<endl;
                        cout<<"j1: "<<(j1pj2[0]+j1mj2[0])/2<<" "<<(j1pj2[1]+j1mj2[1])/2<<" "<<(j1pj2[2]+j1mj2[2])/2<<" "<<(j1pj2[3]+j1mj2[3])/2<<" "<<endl;
                        cout<<"j2: "<<(j1pj2[0]-j1mj2[0])/2<<" "<<(j1pj2[1]-j1mj2[1])/2<<" "<<(j1pj2[2]-j1mj2[2])/2<<" "<<(j1pj2[3]-j1mj2[3])/2<<" "<<endl;
                    };

                dreal hs=1;
                vector<vector<dreal>> jacques_vecs={j1pj2,j1pj2,j1pj2,j1pj2,j1mj2,j1mj2,j1mj2,final_parton};

                dreal tLz=b_flow.t_val(jacques_vecs[0],3);
                jacques_vecs={b_flow.perform_flow(jacques_vecs[0],3,tLz),b_flow.perform_flow(jacques_vecs[1],3,tLz),b_flow.perform_flow(jacques_vecs[2],3,tLz),b_flow.perform_jacques(jacques_vecs[3],tLz,3),b_flow.perform_flow(jacques_vecs[4],3,tLz),b_flow.perform_flow(jacques_vecs[5],3,tLz),b_flow.perform_flow(jacques_vecs[6],3,tLz),b_flow.perform_flow(jacques_vecs[7],3,tLz)};


                dreal tLy=b_flow.t_val(jacques_vecs[0],2);
                jacques_vecs={b_flow.perform_flow(jacques_vecs[0],2,tLy),b_flow.perform_flow(jacques_vecs[1],2,tLy),b_flow.perform_jacques(jacques_vecs[2],tLy,2),b_flow.perform_flow(jacques_vecs[3],2,tLy),b_flow.perform_flow(jacques_vecs[4],2,tLy),b_flow.perform_flow(jacques_vecs[5],2,tLy),b_flow.perform_flow(jacques_vecs[6],2,tLy),b_flow.perform_flow(jacques_vecs[7],2,tLy)};


                dreal tLx=b_flow.t_val(jacques_vecs[0],1);
                jacques_vecs={b_flow.perform_flow(jacques_vecs[0],1,tLx),b_flow.perform_jacques(jacques_vecs[1],tLx,1),b_flow.perform_flow(jacques_vecs[2],1,tLx),b_flow.perform_flow(jacques_vecs[3],1,tLx),b_flow.perform_flow(jacques_vecs[4],1,tLx),b_flow.perform_flow(jacques_vecs[5],1,tLx),b_flow.perform_flow(jacques_vecs[6],1,tLx),b_flow.perform_flow(jacques_vecs[7],1,tLx)};


                dreal txy=r_flow.t_val(jacques_vecs[4],1);
                jacques_vecs={r_flow.perform_flow(jacques_vecs[0],1,txy),r_flow.perform_flow(jacques_vecs[1],1,txy),r_flow.perform_flow(jacques_vecs[2],1,txy),r_flow.perform_flow(jacques_vecs[3],1,txy),r_flow.perform_flow(jacques_vecs[4],1,txy),r_flow.perform_jacques(jacques_vecs[5],txy,1),r_flow.perform_flow(jacques_vecs[6],1,txy),r_flow.perform_flow(jacques_vecs[7],1,txy)};


                dreal tyz=r_flow.t_val(jacques_vecs[4],2);
                jacques_vecs={r_flow.perform_flow(jacques_vecs[0],2,tyz),r_flow.perform_flow(jacques_vecs[1],2,tyz),r_flow.perform_flow(jacques_vecs[2],2,tyz),r_flow.perform_flow(jacques_vecs[3],2,tyz),r_flow.perform_flow(jacques_vecs[4],2,tyz),r_flow.perform_flow(jacques_vecs[5],2,tyz),r_flow.perform_jacques(jacques_vecs[6],tyz,2),r_flow.perform_flow(jacques_vecs[7],2,tyz)};


                dreal tLz_gen=-a/sqrt(1+pow(a,2.));
                jacques_vecs={b_flow.perform_flow(jacques_vecs[0],3,tLz_gen),b_flow.perform_flow(jacques_vecs[1],3,tLz_gen),b_flow.perform_flow(jacques_vecs[2],3,tLz_gen),b_flow.perform_flow(jacques_vecs[3],3,tLz_gen),b_flow.perform_flow(jacques_vecs[4],3,tLz_gen),b_flow.perform_flow(jacques_vecs[5],3,tLz_gen),b_flow.perform_flow(jacques_vecs[6],3,tLz_gen),b_flow.perform_flow(jacques_vecs[7],3,tLz_gen)};


                /*jacques_vecs={b_flow.perform_flow(jacques_vecs[0],3,tLz),b_flow.perform_flow(jacques_vecs[1],3,tLz),b_flow.perform_flow(jacques_vecs[2],3,tLz),b_flow.perform_flow(jacques_vecs[3],3,tLz),b_flow.perform_flow(jacques_vecs[4],3,tLz),b_flow.perform_flow(jacques_vecs[5],3,tLz)};*/


                dreal res_f=Sqrt(SP(jacques_vecs[0],jacques_vecs[0]))*(b_flow.h(tLx)*b_flow.h(tLy)*b_flow.h(tLz)*r_flow.h(txy)*r_flow.h(tyz))/abs(jacques_vecs[1][1]*jacques_vecs[2][2]*jacques_vecs[3][3]*jacques_vecs[5][1]*jacques_vecs[6][2]);
                

                if(DEBUG==4){
                    for(int i=0; i<7; i++){
                        cout<<"vec "<<i<< ":   ";
                        for(int j=0; j<4; j++){
                            cout<<jacques_vecs[i][j]<<" ";
                        };
                        cout<<endl;
                    };

                    cout<<"--------------------"<<endl;

                    cout<<"ts"<<endl;
                    cout<<"Lx: "<<tLx<<endl;
                    cout<<"Ly: "<<tLy<<endl;
                    cout<<"Lz: "<<tLz<<endl;
                    cout<<"Rxy: "<<txy<<endl;
                    cout<<"Ryz: "<<tyz<<endl;

                    cout<<"jacques transf"<<endl;
                    cout<<"Lx: "<<jacques_vecs[1][1]<<endl;
                    cout<<"Ly: "<<jacques_vecs[2][2]<<endl;
                    cout<<"Lz: "<<jacques_vecs[3][3]<<endl;
                    cout<<"Rxy: "<<jacques_vecs[5][1]<<endl;
                    cout<<"Ryz: "<<jacques_vecs[6][2]<<endl;

                    cout<<"hts"<<endl;
                    cout<<"Lx: "<<b_flow.h(tLx)<<endl;
                    cout<<"Ly: "<<b_flow.h(tLy)<<endl;
                    cout<<"Lz: "<<b_flow.h(tLz)<<endl;
                    cout<<"Rxy: "<<r_flow.h(txy)<<endl;
                    cout<<"Ryz: "<<r_flow.h(tyz)<<endl;

                    cout<<"jac   "<<res_f<<endl;
                };

                //dreal x1px2=jacques_vecs[0][0]/eCM;
                //dreal x1mx2=jacques_vecs[3][0]/eCM;

                dreal x1px2=2*(jacques_vecs[0][0])/eCM;
                dreal x1mx2=2*(jacques_vecs[0][3])/eCM;
                dreal pt=Sqrt(pow(jacques_vecs[7][1],2.)+pow(jacques_vecs[7][2],2.));

                if(DEBUG==1){
                    cout<<"j1+j2 post: "<<jacques_vecs[0][0]<<" "<<jacques_vecs[0][1]<<" "<<jacques_vecs[0][2]<<" "<<jacques_vecs[0][3]<<" "<<endl;
                    cout<<"j1-j2 post: "<<jacques_vecs[4][0]<<" "<<jacques_vecs[4][1]<<" "<<jacques_vecs[4][2]<<" "<<jacques_vecs[4][3]<<" "<<endl;
                    cout<<"l post: "<<jacques_vecs[7][0]<<" "<<jacques_vecs[7][1]<<" "<<jacques_vecs[7][2]<<" "<<jacques_vecs[7][3]<<" "<<endl;
                    cout<<"j1 post: "<<(jacques_vecs[0][0]+jacques_vecs[4][0])/2<<" "<<(jacques_vecs[0][1]+jacques_vecs[4][1])/2<<" "<<(jacques_vecs[0][2]+jacques_vecs[4][2])/2<<" "<<(jacques_vecs[0][3]+jacques_vecs[4][3])/2<<" "<<endl;
                    cout<<"j2 post: "<<(jacques_vecs[0][0]-jacques_vecs[4][0])/2<<" "<<(jacques_vecs[0][1]-jacques_vecs[4][1])/2<<" "<<(jacques_vecs[0][2]-jacques_vecs[4][2])/2<<" "<<(jacques_vecs[0][3]-jacques_vecs[4][3])/2<<" "<<endl;
                    cout<<"x1: "<<(x1mx2+x1px2)/2<<" x2: "<<(x1px2-x1mx2)/2<<endl;
                };

                observable res_vec;

                vector<dreal> j1={(jacques_vecs[0][0]+jacques_vecs[4][0])/2,(jacques_vecs[0][1]+jacques_vecs[4][1])/2,(jacques_vecs[0][2]+jacques_vecs[4][2])/2,(jacques_vecs[0][3]+jacques_vecs[4][3])/2};
                vector<dreal> j2={(jacques_vecs[0][0]-jacques_vecs[4][0])/2,(jacques_vecs[0][1]-jacques_vecs[4][1])/2,(jacques_vecs[0][2]-jacques_vecs[4][2])/2,(jacques_vecs[0][3]-jacques_vecs[4][3])/2};
                vector<dreal> pg={jacques_vecs[7][0],jacques_vecs[7][1],jacques_vecs[7][2],jacques_vecs[7][3]};

                res_vec.eval=0;
                res_vec.jac=res_f;
                /*res_vec.j1=j1;
                res_vec.j2=j2;*/
                res_vec.j1=jacques_vecs[0];
                res_vec.j2=jacques_vecs[4];
                res_vec.spin1=double(tot_spin1);
                res_vec.spin2=double(tot_spin2);
                res_vec.pg=pg;

                if(DEBUG==1){
                    cout<<"in clustering"<<endl;
                    cout<<res_vec.eval<<endl;
                    cout<<"jacques: "<<res_vec.jac<<endl;
                    cout<<j1[0]<<" "<<j1[1]<<" "<<j1[2]<<" "<<j1[3]<<endl;
                    cout<<j2[0]<<" "<<j2[1]<<" "<<j2[2]<<" "<<j2[3]<<endl;
                };

                return res_vec;

            };

            vector<dreal> j1={0,0,0,0};
            vector<dreal> j2={0,0,0,0};
            vector<dreal> pg={0,0,0,0};

            observable res_vec_n;

            res_vec_n.eval=0;
            res_vec_n.jac=0;
            res_vec_n.j1=j1;
            res_vec_n.j2=j2;
            res_vec_n.spin1=0;
            res_vec_n.spin2=0;
            res_vec_n.pg=pg;
            
            return res_vec_n;

        };


};




class deformation_field{

    private:

        dreal lambda;
        dreal mij;
        dreal branch_cut_lambda;

    public:

        deformation_field(dreal nlambda, dreal nmij, dreal nbranch_cut_lambda){
            lambda=nlambda;
            mij=nmij;
            branch_cut_lambda=nbranch_cut_lambda;
        };

        void set_hyperparameters(dreal nlambda, dreal nmij, dreal nbranch_cut_lambda){
            lambda=nlambda;
            mij=nmij;
            branch_cut_lambda=nbranch_cut_lambda;
        };


        /*vector<dreal> get_deformation_dt(vector<dreal> k, vector<dreal> l, vector<dreal> q){

            dreal e_k=norm(k);
            dreal e_kl=norm(vector_plus(k,vector_swap(l)));
            dreal e_l=norm(l);
            dreal e_kq=norm(vector_plus(k,vector_swap(q)));
            dreal e_lq=norm(vector_plus(l,vector_swap(q)));

            dreal pinch1=e_k+e_kl-e_l;
            dreal pinch2=e_kq+e_kl-e_lq;

            dreal tp1=pow(pinch1,2.)/(pow(pinch1,2.)+pow(mij,2.));
            dreal tp2=pow(pinch2,2.)/(pow(pinch2,2.)+pow(mij,2.));

            vector<dreal> defo=vector_times(k, tp1*tp2);
            dreal e_defo=norm(defo);

            if(e_kl/e_defo<e_kq/e_defo && e_kl/e_defo<branch_cut_lambda){
                return vector_times(k,lambda*e_kl/e_k);
            };
            if(e_kq/e_defo<e_kl/e_defo && e_kq/e_defo<branch_cut_lambda){
                return vector_times(k,e_kq/e_k*lambda);
            };
            if(branch_cut_lambda<e_kq/e_defo && branch_cut_lambda<e_kl/e_defo){
                return vector_times(defo, branch_cut_lambda*lambda);
            };

            return k;

        };*/


        dual<dreal> get_deformation_dt(vector<dreal> k, vector<dreal> l, vector<dreal> q, int i, int j){

            dual<dreal> kx=k[0]; dual<dreal> ky=k[1]; dual<dreal> kz=k[2];
            dreal lx=l[0]; dreal ly=l[1]; dreal lz=l[2];
            dreal qx=q[0]; dreal qy=q[1]; dreal qz=q[2]; 

            if(i==1){
                kx=k[0]+1_e;
            };
            if(i==2){
                ky=k[1]+1_e; 
            };
            if(i==3){
                kz=k[2]+1_e;
            };


            dual<dreal> e_k=sqrt(kx*kx+ky*ky+kz*kz);
            dual<dreal> e_kl=sqrt((kx-lx)*(kx-lx)+(ky-ly)*(ky-ly)+(kz-lz)*(kz-lz));
            dreal e_l=norm(l);
            dual<dreal> e_kq=sqrt((kx-qx)*(kx-qx)+(ky-qy)*(ky-qy)+(kz-qz)*(kz-qz));
            dreal e_lq=norm(vector_plus(l,vector_swap(q)));

            dual<dreal> pinch1=e_k+e_kl-e_l;
            dual<dreal> pinch2=e_kq+e_kl-e_lq;

            dual<dreal> tp1=pow(pinch1,2.)/(pow(pinch1,2.)+pow(mij,2.));
            dual<dreal> tp2=pow(pinch2,2.)/(pow(pinch2,2.)+pow(mij,2.));

            if(e_kl/e_k<e_kq/e_k && e_kl/e_k<branch_cut_lambda){
                if(j==1){return k[0]*tp1*tp2*e_kl/e_k;};
                if(j==2){return k[1]*tp1*tp2*e_kl/e_k;};
                if(j==3){return k[2]*tp1*tp2*e_kl/e_k;};
            }
            if(e_kq/e_k<e_kl/e_k && e_kq/e_k<branch_cut_lambda){
                if(j==1){return k[0]*tp1*tp2*e_kq/e_k;};
                if(j==2){return k[1]*tp1*tp2*e_kq/e_k;};
                if(j==3){return k[2]*tp1*tp2*e_kq/e_k;};
            }
            if(branch_cut_lambda<e_kq/e_k && branch_cut_lambda<e_kl/e_k){
                if(j==1){return k[0]*tp1*tp2*branch_cut_lambda;};
                if(j==2){return k[1]*tp1*tp2*branch_cut_lambda;};
                if(j==3){return k[2]*tp1*tp2*branch_cut_lambda;};
            }

            dreal c{0};

            return c;

        };


        dcomp jacques_deformation_dt(vector<dreal> k, vector<dreal> l, vector<dreal> q, int signn){

            dreal sign=signn;

            dcomp i{0,1};
            dcomp a11=1.+sign*i*get_deformation_dt(k,l,q,1,1).dpart(); dcomp a12=sign*i*get_deformation_dt(k,l,q,1,2).dpart(); dcomp a13=sign*i*get_deformation_dt(k,l,q,1,3).dpart();
            dcomp a21=sign*i*get_deformation_dt(k,l,q,2,1).dpart(); dcomp a22=1.+sign*i*get_deformation_dt(k,l,q,2,2).dpart(); dcomp a23=sign*i*get_deformation_dt(k,l,q,2,3).dpart();
            dcomp a31=sign*i*get_deformation_dt(k,l,q,3,1).dpart(); dcomp a32=sign*i*get_deformation_dt(k,l,q,3,2).dpart(); dcomp a33=1.+sign*i*get_deformation_dt(k,l,q,3,3).dpart();

            dcomp jacques=a11*(a22*a33-a23*a32)-a12*(a21*a33-a23*a31)+a13*(a21*a32-a22*a31);

            return jacques;

        };


};



dreal group_average(int spin1, int spin2){
    int spins[2]={spin1, spin2};
    dreal res=1;
    for(int i=0;i<2;i++){
        if(spins[i]==0){
            res*=1/8.;
        };
        if(abs(spins[i])>0){
            res*=1/(pow(3.,abs(spins[i])));
        };
    };
    return res;
};


class integrands{

    private:
        
        dreal m;
        dreal eCM;
        dreal constant;
        dreal MUV;

        observables my_obs;
        causal_flow my_flow;
        deformation_field my_defo;

    public:
        integrands(dreal nm, dreal neCM, dreal nsigma, dreal nMUV):my_obs(0,0,neCM,nm), my_flow(nm,nsigma), my_defo(-1,0,-1){

            
            eCM=neCM;
            m=nm;
            MUV=nMUV;
            dreal gV=0.1303037631031056;
            dreal gs=sqrt(0.118*4*M_PI);
            int NC=3;
            dreal CF=(NC*NC-1)/(2*NC);
            dreal spin_norm=2;

            dreal PS_constant=1/(16*pow(M_PI,2.)*eCM*eCM);

            constant=pow(gV,2)*pow(gs,2)*CF*NC/(pow(spin_norm,2.))*PS_constant*0.389379*pow(10,9);
            
    
            
        }

        integrands(dreal nm, dreal neCM, dreal nsigma, dreal nMUV, dreal res_c, dreal res_s):my_obs(res_c,res_s,neCM,nm), my_flow(nm,nsigma), my_defo(-1,0,-1){
            eCM=neCM;
            m=nm;
            /*dreal gV=0.1303037631031056;
            dreal gs=Sqrt(0.118*4*M_PI);
            int NC=3;
            dreal CF=(NC*NC-1)/(2*NC);
            dreal spin_norm=2;

            dreal PS_constant=1/(16*pow(M_PI,2.)*eCM*eCM);
            constant=pow(gV,2)*pow(gs,2)*CF*NC/(pow(spin_norm,2.))*PS_constant*0.389379*pow(10,9);
            cout<<constant<<endl;*/

            dreal gV=0.1303037631031056;
            dreal gs=sqrt(0.118*4*M_PI);
            dreal NC=3;
            dreal CF=(NC*NC-1)/(2*NC);
            dreal spin_norm=2;
            dreal consty=pow(gV,2)*pow(gs,2)*CF*NC/(pow(spin_norm,2.)*16*pow(M_PI,2.));


        //dreal constant_NLO=-pow(0.1303037631031056,2)*0.389379*pow(10,9)*(0.118*4*M_PI)*4/(3*3*2*2*4*eCM*eCM)/pow(2*M_PI,3.);
            constant=consty/(eCM*eCM)*0.389379*pow(10,9);


            MUV=nMUV;
        }
        
        void set_MUV(dreal nMUV){MUV=nMUV;};

        void set_kinematics(dreal nm, dreal neCM){my_obs.set_kinematics(nm,neCM); m=nm; eCM=neCM;};

        void set_res(dreal nres_c, dreal nres_s){my_obs.set_res(nres_c,nres_s);};

        void set_sigma(dreal nsigma){my_flow.set_parameters(m, nsigma);};

        void set_deformation_hyperparameters(dreal nlambda, dreal nmij, dreal nbc_lambda){
            my_defo.set_hyperparameters(nlambda, nmij, nbc_lambda);
        };



        dreal LO_scalar(vector<dreal> p, vector<dreal> q, dreal a){
            
            vector<vector<dreal>> ni={p,vector_minus(q,p)};
            
            vector<vector<dreal>> nf={};

            vector<dreal> ps=my_flow.perform_flow(ni,nf,q,p);
            vector<dreal> qs=my_flow.perform_flow(ni,nf,q,q);

            dreal jac_flow=my_flow.jacques(ni,nf,q,2);
            dreal flow_h=my_flow.h(my_flow.t_val(ni,nf,q));

            dreal e1=norm(ps);
            dreal e2=sqrt(pow(norm(qs),2.)+pow(m,2));
            dreal e3=norm(vector_minus(ps,qs));

            ps.insert(ps.begin(), e1);
            qs.insert(qs.begin(), e2);

            if(DEBUG==1){
                cout<<"pmq"<<endl;
                cout<<vector_minus(q,p)[0]<<" "<<vector_minus(q,p)[1]<<" "<<vector_minus(q,p)[2]<<endl;
                cout<<"t"<<endl;
                cout<<my_flow.t_val(ni,nf,q)<<endl;
                cout<<"flow jacobian:  "<<jac_flow<<endl;
                cout<<"flow h function:  "<<flow_h<<endl;
                cout<<"pscaled"<<endl;
                cout<<ps[0]<<" "<<ps[1]<<" "<<ps[2]<<endl;
                cout<<"qscaled"<<endl;
                cout<<qs[0]<<" "<<qs[1]<<" "<<qs[2]<<endl;
                cout<<"p4scaled"<<endl;
                cout<<ps[0]<<" "<<ps[1]<<" "<<ps[2]<<" "<<ps[3]<<endl;
                cout<<"q4scaled"<<endl;
                cout<<qs[0]<<" "<<qs[1]<<" "<<qs[2]<<" "<<qs[3]<<endl;
                cout<<"(q4-p4)scaled"<<endl;
                cout<<vector_minus(qs,ps)[0]<<" "<<vector_minus(qs,ps)[1]<<" "<<vector_minus(qs,ps)[2]<<" "<<vector_minus(qs,ps)[3]<<endl;
            };

            vector<vector<dreal>> constituents={ps, vector_minus(qs,ps)};

            int spins[2]={1,1};

            vector<dreal> obs=my_obs.b2b_sampling(constituents, 2, spins, a);

            return jac_flow*flow_h*obs[0]/(8*e1*e2*e3);

        };

        dreal LO_DrellYan_numerator(vector<dreal> p, vector<dreal> q){
            return 8*(SP(p,p)-SP(p,q));
        };

        dreal LO_DrellYan_integrand(vector<dreal> p, vector<dreal> q, dreal a){
            vector<vector<dreal>> ni={p,vector_minus(q,p)};
            
            vector<vector<dreal>> nf={};

            vector<dreal> ps=my_flow.perform_flow(ni,nf,q,p);
            vector<dreal> qs=my_flow.perform_flow(ni,nf,q,q);

            dreal jac_flow=my_flow.jacques(ni,nf,q,2);
            dreal flow_h=my_flow.h(my_flow.t_val(ni,nf,q));

            dreal e1=norm(ps);
            dreal e2=sqrt(pow(norm(qs),2.)+pow(m,2));
            dreal e3=norm(vector_minus(ps,qs));

            ps.insert(ps.begin(), e1);
            qs.insert(qs.begin(), e2);

            vector<vector<dreal>> constituents={ps, vector_minus(qs,ps)};

            int spins[2]={1,1};

            vector<dreal> obs_res=my_obs.b2b_sampling(constituents, 2, spins, a);;

            if(DEBUG==1){
                    vector<dreal> psp={0.0, 0.0, 91.188/2};
                    vector<dreal> qsp={0.0, 0.0, 0.0};

                    dreal epp=norm(psp);
                    dreal eqp=sqrt(pow(norm(qsp),2.)+pow(m,2));
                    dreal epq=norm(psp);


                    psp.insert(psp.begin(), epp);
                    qsp.insert(qsp.begin(), eqp);

                    cout<<"LO BORN: "<<LO_DrellYan_numerator(psp,qsp)/(8*epp*eqp*epq)<<endl;
            };

            return constant*jac_flow*flow_h*obs_res[0]*LO_DrellYan_numerator(ps,qs)/(8*e1*e2*e3);
        };





        dreal NLO_DT_DrellYan_numerator(vector<dreal> k, vector<dreal> l, vector<dreal> q){
            return 128*(-pow(SP(k, l),2.) + SP(k, l)*SP(k, q) + SP(k, l)*SP(l, q) - SP(k, q)*SP(l, q));
        };

        dcomp NLO_DT_DrellYan_numerator(vector<dcomp> k, vector<dreal> l, vector<dreal> q){
            return 128.*(-pow(SP(k, l),2.) + SP(k, l)*SP(k, q) + SP(k, l)*SP(l, q) - SP(k, q)*SP(l, q));
        };


        //1/(MUV^2 + norm[k]^2)^(5/2)
        dreal NLO_DT_UV_DrellYan_numeratorA(vector<dreal> k, vector<dreal> l, vector<dreal> q){
            return -24*pow(SP(k, l),2.) + 24*SP(k, l)*SP(k, q) + 12*pow(MUV,2.)*SP(l, l) + 12*SP(k, k)*SP(l, l) - 
                12*pow(MUV,2.)*SP(l, q) - 12*SP(k, k)*SP(l, q);
        };

        //1/(MUV^2 + norm[k]^2)^(3/2)
        dreal NLO_DT_UV_DrellYan_numeratorB(vector<dreal> k, vector<dreal> l, vector<dreal> q){
            return 12*SP(l, l) - 12*SP(l, q) + 8*pow(l[0],2) - 8*l[0]*q[0]; 
        };

        dcomp NLO_DT_UV_DrellYan_numeratorA(vector<dcomp> k, vector<dreal> l, vector<dreal> q){
            return -24.*pow(SP(k, l),2.) + 24.*SP(k, l)*SP(k, q) + 12.*pow(MUV,2.)*SP(l, l) + 12.*SP(k, k)*SP(l, l) - 
                12*pow(MUV,2.)*SP(l, q) - 12.*SP(k, k)*SP(l, q);
        };

        //1/(MUV^2 + norm[k]^2)^(3/2)
        dcomp NLO_DT_UV_DrellYan_numeratorB(vector<dcomp> k, vector<dreal> l, vector<dreal> q){
            return 12.*SP(l, l) - 12.*SP(l, q) + 8.*pow(l[0],2) - 8.*l[0]*q[0]; 
        };


        dreal NLO_SE_DrellYan_numerator(vector<dreal> k, vector<dreal> l, vector<dreal> q){
            return 64*(SP(k, k)*SP(k, l) - 2*SP(k, l)*SP(k, q) + SP(k, k)*SP(l, q));
        };

        dcomp NLO_SE_DrellYan_numerator(vector<dreal> k, vector<dcomp> l, vector<dreal> q){
            return 64.*(SP(k, k)*SP(k, l) - 2.*SP(k, l)*SP(k, q) + SP(k, k)*SP(l, q));
        };


        //THIS DERIVATIVE LOOKS VERY WRONG (USING ROUTING OF ORIGINAL BUBBLE, BUT DERIVATIVE OF TARGET DIAGRAM
        dreal NLO_SE_DrellYan_numerator_d(vector<dreal> k, vector<dreal> l, vector<dreal> q){
            return 64*(SP(k, k)*k[0] - 2*k[0]*SP(k, q) + SP(k, k)*q[0]);
        };  

        dcomp NLO_SE_DrellYan_numerator_d(vector<dreal> k, vector<dcomp> l, vector<dreal> q){
            return 64*(SP(k, k)*k[0] - 2*k[0]*SP(k, q) + SP(k, k)*q[0]);
        };  


        //ADD ep+em->d+dbar flow

        observable_c NLO_dt_DrellYan_integrand(vector<dreal> k, vector<dreal> l, vector<dreal> q, dreal a, int comp){


            if(comp==1 || comp==2){


                vector<dreal> l_in=l;
                vector<dreal> k_in=k;
                vector<dreal> q_in=q;

                dreal defo_sign;
                dreal jac_sign;

                if(comp==1){
                    l_in=l;
                    k_in=k;
                    q_in=q;
                    defo_sign=1;
                    jac_sign=1;
                }else{
                    l_in=k;//vector_swap(k);
                    k_in=l;//vector_swap(l);
                    q_in=q;//vector_swap(q);
                    defo_sign=-1;
                    jac_sign=-1;
                };

                vector<vector<dreal>> ni={l_in,vector_swap(vector_minus(l_in,q_in))};
                vector<vector<dreal>> nf={};

                if(my_flow.t_val(ni,nf,q_in)==999){
                    observable_c res;
                    vector<dreal> j1={0,0,0,0};
                    vector<dreal> j2={0,0,0,0};
                    vector<dreal> pg={0,0,0,0};
                    res.eval=0;
                    res.jac=0;
                    res.j1=j1;
                    res.j2=j2;
                    res.spin1=0;
                    res.spin2=0;
                    res.pg=pg;
                    return res;
                };
                
                vector<dreal> ks_real=my_flow.perform_flow(ni,nf,q_in,k_in);
                vector<dreal> ls=my_flow.perform_flow(ni,nf,q_in,l_in);
                vector<dreal> qs=my_flow.perform_flow(ni,nf,q_in,q_in);

                dreal jac_flow=my_flow.jacques(ni,nf,q_in,3);
                dreal flow_h=my_flow.h(my_flow.t_val(ni,nf,q_in));

                dcomp i(0,1.);
                //vector<dreal> kappa=my_defo.get_deformation(ks_real,ls,1);
                vector<dreal> kappa={my_defo.get_deformation_dt(ks_real,ls,qs,1,1).rpart(),my_defo.get_deformation_dt(ks_real,ls,qs,1,2).rpart(),my_defo.get_deformation_dt(ks_real,ls,qs,1,3).rpart()};
                vector<dcomp> ks={ks_real[0]+defo_sign*i*kappa[0],ks_real[1]+defo_sign*i*kappa[1],ks_real[2]+defo_sign*i*kappa[2]};
                //dcomp defo_jacques=my_defo.jacques_deformation(ks_real,ls,1,jac_sign);
                dcomp defo_jacques=my_defo.jacques_deformation_dt(ks_real,ls,qs,double(jac_sign));

                /*cout<<"k: "<<ks_real[0]<<" "<<ks_real[1]<<" "<<ks_real[2]<<endl;
                cout<<"l: "<<ls[0]<<" "<<ls[1]<<" "<<ls[2]<<endl;
                cout<<"q: "<<qs[0]<<" "<<qs[1]<<" "<<qs[2]<<endl;
                cout<<"kappa:  "<<kappa[0]<<" "<<kappa[1]<<" "<<kappa[2]<<endl;
                cout<<"jacques:  "<<defo_jacques<<endl;*/

                dcomp e_k=norm(ks);
                dcomp e_kq=norm(vector_minus(ks,qs));
                dreal e_l=norm(ls);
                dcomp e_kl=norm(vector_minus(ks,ls));
                dreal e_q=sqrt(pow(norm(qs),2.)+pow(m,2));

                dcomp e_kUV=Sqrt(pow(norm(ks),2.)+pow(MUV,2));
                
                
                ls.insert(ls.begin(), e_l);
                qs.insert(qs.begin(), e_q);

                vector<dcomp> ks1={e_k,ks[0],ks[1],ks[2]};

                dcomp dt_virtual_left_1=NLO_DT_DrellYan_numerator(ks1,ls,qs)/
                    (2.*e_l*2.*e_k*2.*(-e_l+e_q)*2.*e_q*SP(vector_minus(ks1,ls),vector_minus(ks1,ls))*SP(vector_minus(ks1,qs),vector_minus(ks1,qs)));

                vector<dcomp> ks2={e_l+e_kl,ks[0],ks[1],ks[2]};

                dcomp dt_virtual_left_2=NLO_DT_DrellYan_numerator(ks2, ls, qs)/
                    ((2.*e_l*2.*(-e_l+e_q)*2.*(e_kl)*2.*e_q)*SP(ks2,ks2)*SP(vector_minus(ks2,qs),vector_minus(ks2,qs)));

                vector<dcomp> ks3={e_kq+e_q,ks[0],ks[1],ks[2]};

                dcomp dt_virtual_left_3=NLO_DT_DrellYan_numerator(ks3, ls, qs)/
                    ((2.*e_l*2.*(-e_l+e_q)*2.*(e_kq)*2.*e_q)*SP(ks3,ks3)*SP(vector_minus(ks3,ls),vector_minus(ks3,ls)));

                vector<dcomp> ksuv={0,ks[0],ks[1],ks[2]};

                dcomp dt_left_uv=1./(2.*e_l*2.*(-e_l+e_q)*2.*e_q)*(NLO_DT_UV_DrellYan_numeratorA(ksuv, ls, qs)/pow(e_kUV,5.)+
                    NLO_DT_UV_DrellYan_numeratorB(ksuv, ls, qs)/pow(e_kUV,3.)
                    );

                /*if(comp==2){
                    cout<<"scaling"<<endl;
                    cout<<"e_l: "<<e_k<<endl;
                    cout<<real(pow(e_k,3.)*dt_virtual_left_1)<<" "<<real(pow(e_k,3.)*dt_virtual_left_2)<<" "<<real(pow(e_k,3.)*dt_virtual_left_3)<<"     tot orig:"<<real(pow(e_k,3.)*(dt_virtual_left_1+dt_virtual_left_2+dt_virtual_left_3-dt_left_uv))<<endl;
                    cout<<"no scaling"<<endl;
                    cout<<real(dt_virtual_left_1)<<" "<<real(dt_virtual_left_2)<<" "<<real(dt_virtual_left_3)<<"     tot orig:"<<real((dt_virtual_left_1+dt_virtual_left_2+dt_virtual_left_3-dt_left_uv))<<endl;
                };*/

                if(DEBUG==3){
                    cout<<"ks1: "<<ks1[0]<<" "<<ks1[1]<<" "<<ks1[2]<<" "<<ks1[3]<<endl;
                    cout<<"ks2: "<<ks2[0]<<" "<<ks2[1]<<" "<<ks2[2]<<" "<<ks2[3]<<endl;
                    cout<<"ks3: "<<ks3[0]<<" "<<ks3[1]<<" "<<ks3[2]<<" "<<ks3[3]<<endl;
                    cout<<"ls: "<<ls[0]<<" "<<ls[1]<<" "<<ls[2]<<" "<<ls[3]<<endl;
                    cout<<"qs: "<<qs[0]<<" "<<qs[1]<<" "<<qs[2]<<" "<<qs[3]<<endl;
                    cout<<SP(vector_minus(ks1,ls),vector_minus(ks1,ls))<<endl;
                    cout<<SP(vector_minus(ks1,qs),vector_minus(ks1,qs))<<endl;
                    cout<<SP(ks2,ks2)<<endl;
                    cout<<SP(vector_minus(ks2,qs),vector_minus(ks2,qs))<<endl;
                    cout<<SP(ks3,ks3)<<endl;
                    cout<<SP(vector_minus(ks3,ls),vector_minus(ks3,ls))<<endl;
                    cout<<"res: "<<dt_left_uv<<endl;
                };

                vector<vector<dreal>> constituents={ls,vector_swap(vector_minus(ls,qs))};


                int spins[2]={1,-1};

                vector<dreal> pg={0,0,0,0};

                observable obs_res_n=my_obs.b2b_sampling_final(constituents, 2, pg, spins, a);

                observable_c obs_res;
                obs_res.eval=obs_res_n.eval;
                obs_res.jac=obs_res_n.jac;
                obs_res.j1=obs_res_n.j1;
                obs_res.j2=obs_res_n.j2;
                obs_res.spin1=obs_res_n.spin1;
                obs_res.spin2=obs_res_n.spin2;
                obs_res.pg=obs_res_n.pg;


                
                if(DEBUG==1){
                    vector<dreal> ksp={0.3, 0.1, 0.5};
                    vector<dreal> lsp={0.2, 0.7, 0.1};
                    vector<dreal> qsp={0.4, 0.4, 0.1};

                    dreal ekp=norm(ksp);
                    dreal ekqp=norm(vector_minus(ksp,qsp));
                    dreal eqp=sqrt(pow(norm(qsp),2.)+pow(m,2));
                    dreal elp=norm(lsp);  
                    dreal eklp=norm(vector_minus(ksp,lsp));  
                    dreal eklqp=norm(vector_plus(lsp,vector_plus(ksp,qsp)));
                    dreal e_kUVp=sqrt(pow(norm(ksp),2.)+pow(MUV,2));

                    lsp.insert(lsp.begin(), elp);
                    qsp.insert(qsp.begin(), eqp);
                    

                    vector<dreal> ks1p={ekp,ksp[0],ksp[1],ksp[2]};

                    dreal dt_virtual_left_1p=NLO_DT_DrellYan_numerator(ks1p,lsp,qsp)/
                        (2.*elp*2.*ekp*2.*(-elp+eqp)*2*eqp*SP(vector_minus(ks1p,lsp),vector_minus(ks1p,lsp))*SP(vector_minus(ks1p,qsp),vector_minus(ks1p,qsp)));
                    

                    vector<dreal> ks2p={elp+eklp,ksp[0],ksp[1],ksp[2]};

                    dreal dt_virtual_left_2p=NLO_DT_DrellYan_numerator(ks2p, lsp, qsp)/
                        ((2.*elp*2.*(-elp+eqp)*2.*(eklp)*2*eqp)*SP(ks2p,ks2p)*SP(vector_minus(ks2p,qsp),vector_minus(ks2p,qsp)));

                    vector<dreal> ks3p={ekqp+eqp,ksp[0],ksp[1],ksp[2]};

                    dreal dt_virtual_left_3p=NLO_DT_DrellYan_numerator(ks3p, lsp, qsp)/
                        ((2.*elp*2.*(-elp+eqp)*2.*(ekqp)*2*eqp)*SP(ks3p,ks3p)*SP(vector_minus(ks3p,lsp),vector_minus(ks3p,lsp)));

                    vector<dreal> ksuvp={0,ksp[0],ksp[1],ksp[2]};

                    dreal dt_left_uv=1/(2.*elp*2.*(-elp+eqp)*2*eqp)*(NLO_DT_UV_DrellYan_numeratorA(ksuvp, lsp, qsp)/pow(e_kUVp,5.)+
                        NLO_DT_UV_DrellYan_numeratorB(ksuvp, lsp, qsp)/pow(e_kUVp,3.)
                        );

                    cout<<NLO_DT_DrellYan_numerator(ks1p,lsp,qsp)<<endl;
                    cout<<"ltd1: "<<dt_virtual_left_1p<<endl;
                    cout<<"energies1: "<<2.*elp*2.*ekp*2.*(-elp+eqp)*2*eqp<<endl;
                    cout<<"sp11: "<<SP(vector_minus(ks1p,lsp),vector_minus(ks1p,lsp))<<endl;
                    cout<<"sp12: "<<SP(vector_minus(ks1p,qsp),vector_minus(ks1p,qsp))<<endl;
                    cout<<"ltd2: "<<dt_virtual_left_2p<<endl;
                    cout<<"ltd3: "<<dt_virtual_left_3p<<endl;
                    cout<<"res: "<<dt_virtual_left_1p+dt_virtual_left_2p+dt_virtual_left_3p<<endl;
                    cout<<"res UV: "<<dt_left_uv<<endl;    

                };

                
                dreal colour_average=group_average(obs_res.spin1,obs_res.spin2);

                //dcomp integrand_val=constant*colour_average*(dt_virtual_left_1+dt_virtual_left_2+dt_virtual_left_3-dt_left_uv);
                dcomp integrand_val=constant*defo_jacques*colour_average*(dt_virtual_left_1+dt_virtual_left_2+dt_virtual_left_3-dt_left_uv);

                //obs_res.jac=obs_res_n.jac*jac_flow*flow_h*defo_jacques;
                obs_res.jac=obs_res_n.jac*jac_flow*flow_h;

                obs_res.eval=integrand_val;

                if(comp==5){
                    cout<<"--------------------"<<endl;
                    cout<<"f64"<<endl;
                    cout<<"eCM "<<dreal(eCM)<<"     mZ: "<<dreal(m)<<endl;
                    cout<<"constant "<<dreal(constant)<<"     MUV: "<<dreal(MUV)<<endl;
                    cout<<"virtual contributions"<<endl;
                    cout<<my_flow.t_val(ni,nf,q_in)<<endl;
                    cout<<jac_flow<<endl;
                    cout<<flow_h<<endl;
                    cout<<obs_res_n.jac<<endl;
                    cout<<defo_jacques<<endl;
                    cout<<"info on virtual2"<<endl;
                    cout<<obs_res.eval<<endl;
                    cout<<obs_res.jac<<endl;
                    cout<<defo_jacques<<endl;
                    cout<<obs_res.j1[0]<<" "<<obs_res.j1[1]<<" "<<obs_res.j1[2]<<" "<<obs_res.j1[3]<<endl;
                    cout<<obs_res.j2[0]<<" "<<obs_res.j2[1]<<" "<<obs_res.j2[2]<<" "<<obs_res.j2[3]<<endl;
                    cout<<obs_res.pg[0]<<" "<<obs_res.pg[1]<<" "<<obs_res.pg[2]<<" "<<obs_res.pg[3]<<endl;
                };


                if(DEBUG==2){
                    cout<<"virtual contributions"<<endl;
                    cout<<my_flow.t_val(ni,nf,q_in)<<endl;
                    cout<<jac_flow<<endl;
                    cout<<flow_h<<endl;
                    cout<<obs_res_n.jac<<endl;
                    cout<<defo_jacques<<endl;
                    cout<<"info on virtual2"<<endl;
                    cout<<obs_res.eval<<endl;
                    cout<<obs_res.jac<<endl;
                    cout<<defo_jacques<<endl;
                    cout<<obs_res.j1[0]<<""<<obs_res.j1[1]<<""<<obs_res.j1[2]<<""<<obs_res.j1[3]<<endl;
                    cout<<obs_res.j2[0]<<""<<obs_res.j2[1]<<""<<obs_res.j2[2]<<""<<obs_res.j2[3]<<endl;
                    cout<<obs_res.pg[0]<<""<<obs_res.pg[1]<<""<<obs_res.pg[2]<<""<<obs_res.pg[3]<<endl;
                };

                return obs_res;

            };



            if(comp==3 || comp==4){

                vector<dreal> l_in=l;
                vector<dreal> k_in=k;
                vector<dreal> q_in=q;

                if(comp==3){
                    l_in=l;
                    k_in=k;
                    q_in=q;
                }else{
                    l_in=vector_swap(k);
                    k_in=vector_swap(l);
                    q_in=vector_swap(q);
                };

                vector<vector<dreal>> ni={k_in,vector_swap(vector_minus(k_in,l_in)),vector_swap(vector_minus(l_in,q_in))};
                vector<vector<dreal>> nf={};

                if(my_flow.t_val(ni,nf,q_in)==999){
                    observable_c res;
                    vector<dreal> j1={0,0,0,0};
                    vector<dreal> j2={0,0,0,0};
                    vector<dreal> pg={0,0,0,0};
                    res.eval=0;
                    res.jac=0;
                    res.j1=j1;
                    res.j2=j2;
                    res.spin1=0;
                    res.spin2=0;
                    res.pg=pg;
                    return res;
                };

                vector<dreal> ks=my_flow.perform_flow(ni,nf,q_in,k_in);
                vector<dreal> ls=my_flow.perform_flow(ni,nf,q_in,l_in);
                vector<dreal> qs=my_flow.perform_flow(ni,nf,q_in,q_in);

                dreal jac_flow=my_flow.jacques(ni,nf,q_in,3);
                dreal flow_h=my_flow.h(my_flow.t_val(ni,nf,q_in));

                dreal e_k=norm(ks);
                dreal e_lq=norm(vector_minus(ls,qs));
                dreal e_kl=norm(vector_minus(ls,ks));
                dreal e_q=sqrt(pow(norm(qs),2.)+pow(m,2));
                
                ks.insert(ks.begin(), e_k);
                ls.insert(ls.begin(), -e_lq+e_q);
                qs.insert(qs.begin(), e_q);

                vector<dreal> ks1={e_k,ks[0],ks[1],ks[2]};

                vector<vector<dreal>> constituents={ks,vector_swap(vector_minus(ks,ls)),vector_swap(vector_minus(ls,qs))};

                int spins[3]={1,0,-1};

                vector<dreal> pg={0,0,0,0};

                observable obs_res_n=my_obs.b2b_sampling_final(constituents, 3, pg, spins, a);

                observable_c obs_res;
                obs_res.eval=obs_res_n.eval;
                obs_res.jac=obs_res_n.jac;
                obs_res.j1=obs_res_n.j1;
                obs_res.j2=obs_res_n.j2;
                obs_res.spin1=obs_res_n.spin1;
                obs_res.spin2=obs_res_n.spin2;
                obs_res.pg=obs_res_n.pg;

                if(DEBUG==1){
                    vector<dreal> ksp={0.3, 0.1, 0.5};
                    vector<dreal> lsp={0.2, 0.7, 0.1};
                    vector<dreal> qsp={0.1, 0.4, 0.5};

                    dreal ekp=norm(ksp);
                    dreal eqp=sqrt(pow(norm(qsp),2.)+pow(m,2));
                    dreal elp=norm(lsp);  
                    dreal ekqp=norm(vector_plus(ksp,qsp));
                    dreal elqp=norm(vector_minus(lsp,qsp)); //HERE  
                    dreal eklqp=norm(vector_plus(lsp,vector_plus(ksp,qsp)));

                    ksp.insert(ksp.begin(), ekp);
                    lsp.insert(lsp.begin(), -elqp+eqp);
                    qsp.insert(qsp.begin(), eqp);


                    cout<<"numerator: "<<NLO_DT_DrellYan_numerator(ksp,lsp,qsp)<<endl;
                    cout<<"energies: "<<2.*ekp*2.*(eqp-elqp-ekp)*2.*(elqp)*2*eqp<<endl;
                    cout<<"sps: "<<SP(vector_minus(ksp,qsp),vector_minus(ksp,qsp))*SP(lsp,lsp)<<endl;
                    cout<<NLO_DT_DrellYan_numerator(ksp,lsp,qsp)/
                    ((2.*ekp*2.*(eqp-elqp-ekp)*2.*(elp)*2*eqp)*
                    SP(vector_minus(ksp,qsp),vector_minus(ksp,qsp))*SP(lsp,lsp))<<endl; 

                };   

                dreal colour_average=group_average(obs_res.spin1,obs_res.spin2);
                    
                dreal integrand_val=constant*colour_average*(NLO_DT_DrellYan_numerator(ks,ls,qs)/
                    ((2.*e_k*2.*(e_q-e_lq-e_k)*2.*(e_lq)*2*e_q)*
                    SP(vector_minus(ks,qs),vector_minus(ks,qs))*SP(ls,ls)
                    ));

                obs_res.jac=obs_res.jac*jac_flow*flow_h;

                obs_res.eval=integrand_val;

                /*cout<<"------------"<<endl;
                cout<<obs_res.eval<<endl;
                cout<<obs_res.jac<<endl;
                cout<<obs_res.spin1<<endl;
                cout<<obs_res_n.j1[0]<<" "<<obs_res_n.j1[1]<<" "<<obs_res_n.j1[2]<<" "<<obs_res_n.j1[3]<<endl;
                cout<<obs_res.spin2<<endl;
                cout<<obs_res_n.j2[0]<<" "<<obs_res_n.j2[1]<<" "<<obs_res_n.j2[2]<<" "<<obs_res_n.j2[3]<<endl;
                cout<<obs_res_n.pg[0]<<" "<<obs_res_n.pg[1]<<" "<<obs_res_n.pg[2]<<" "<<obs_res_n.pg[3]<<endl;*/

                return obs_res;

            };

            observable_c obs_res_f;
            obs_res_f.eval=0;
            obs_res_f.jac=0;
            obs_res_f.j1={0,0,0,0};
            obs_res_f.j2={0,0,0,0};
            obs_res_f.spin1=0;
            obs_res_f.spin2=0;
            obs_res_f.pg={0,0,0,0};


            return obs_res_f;


        };



        vector<dreal> NLO_dt_jj_integrand(vector<dreal> k, vector<dreal> l, dreal q0, int comp, int phase){

            if(comp==1 || comp==2){

                vector<dreal> l_in=l;
                vector<dreal> k_in=k;
                vector<dreal> q_in={q0,0,0,0};

                dreal defo_sign;

                if(comp==1){
                    l_in=l;
                    k_in=k;
                    defo_sign=1;
                }else{
                    l_in=vector_swap(k);
                    k_in=vector_swap(l);
                    //q_in=vector_swap(q_in);
                    defo_sign=-1;
                };

                vector<vector<dreal>> ni={l_in,vector_swap(l_in)};

                if(DEBUG==1){
                    cout<<"----"<<endl;
                    cout<<"presum: "<<2*norm(l_in)<<endl;
                };

                vector<dreal> ks=my_flow.perform_flow_jj(ni,q0,k_in);
                vector<dreal> ls=my_flow.perform_flow_jj(ni,q0,l_in);

                if(DEBUG==1){
                    cout<<ks[0]<<" "<<ks[1]<<" "<<ks[2]<<endl;
                    cout<<ls[0]<<" "<<ls[1]<<" "<<ls[2]<<endl;
                    cout<<ks[0]-ls[0]<<" "<<ks[1]-ls[1]<<" "<<ks[2]-ls[2]<<endl;
                };

                dreal jac_flow=my_flow.jacques_jj(ni,q0,2);
                dreal flow_h=my_flow.h(my_flow.t_val_jj(ni,q0));


                /*cout<<"t value:  "<<my_flow.t_val_jj(ni,q0)<<endl;
                cout<<"t jac value:  "<<my_flow.jacques_jj(ni,q0,2)<<endl;*/

                //dcomp i(0,1.);
                //vector<dreal> kappa=my_defo.get_deformation(ks_real,ls,1);
                //vector<dcomp> ks={ks_real[0]+defo_sign*i*kappa[0],ks_real[1]+defo_sign*i*kappa[1],ks_real[2]+defo_sign*i*kappa[2]};
                //dcomp defo_jacques=my_defo.jacques_deformation(ks_real,ls,1,defo_sign);

                dcomp e_k=norm(ks);
                dcomp e_kq=norm(ks);
                dreal e_l=norm(ls);
                dcomp e_kl=norm(vector_minus(ks,ls));

                dcomp e_kUV=Sqrt(pow(norm(ks),2.)+pow(MUV,2));
                
                
                ls.insert(ls.begin(), e_l);
                

                vector<dcomp> ks1={e_k,ks[0],ks[1],ks[2]};

                dcomp dt_virtual_left_1=NLO_DT_DrellYan_numerator(ks1,ls,q_in)/
                    (2.*e_l*2.*e_k*2.*(-e_l+q_in[0])*SP(vector_minus(ks1,ls),vector_minus(ks1,ls))*SP(vector_minus(ks1,q_in),vector_minus(ks1,q_in)));

                //cout<<"num1: "<<NLO_DT_DrellYan_numerator(ks1,ls,q_in).real()<<endl;

                vector<dcomp> ks2={e_l+e_kl,ks[0],ks[1],ks[2]};

                dcomp dt_virtual_left_2=NLO_DT_DrellYan_numerator(ks2, ls, q_in)/
                    ((2.*e_l*2.*(-e_l+q_in[0])*2.*(e_kl))*SP(ks2,ks2)*SP(vector_minus(ks2,q_in),vector_minus(ks2,q_in)));

                //cout<<"num2: "<<NLO_DT_DrellYan_numerator(ks2, ls, q_in).real()<<endl;

                vector<dcomp> ks3={e_kq+q_in[0],ks[0],ks[1],ks[2]};

                dcomp dt_virtual_left_3=NLO_DT_DrellYan_numerator(ks3, ls, q_in)/
                    ((2.*e_l*2.*(-e_l+q_in[0])*2.*(e_kq))*SP(ks3,ks3)*SP(vector_minus(ks3,ls),vector_minus(ks3,ls)));

                //cout<<"num3: "<<NLO_DT_DrellYan_numerator(ks3, ls, q_in).real()<<endl;

                vector<dcomp> ksuv={0,ks[0],ks[1],ks[2]};

                dcomp dt_left_uv=1./(2.*e_l*2.*(-e_l+q_in[0]))*(NLO_DT_UV_DrellYan_numeratorA(ksuv, ls, q_in)/pow(e_kUV,5.)+
                    NLO_DT_UV_DrellYan_numeratorB(ksuv, ls, q_in)/pow(e_kUV,3.)
                    );

                dreal integrand_val;


                //cout<<"res: "<<dt_virtual_left_1.real()+dt_virtual_left_2.real()+dt_virtual_left_3.real()<<endl;

                dcomp rr=dt_virtual_left_1+dt_virtual_left_2+dt_virtual_left_3+dt_left_uv;

                if(DEBUG==1){
                    cout<<"############VIRTUAL CONTRIBUTIONS########"<<endl;
                    cout<<"q0: "<<q0<<endl;
                    cout<<2*e_l<<endl;
                    cout<<"t: "<<my_flow.t_val_jj(ni,q0)<<endl;
                    cout<<"num1:  "<<NLO_DT_DrellYan_numerator(ks1,ls,q_in)<<endl;
                    dcomp den_en1=2.*e_l*2.*e_k*2.*(-e_l+q_in[0]);
                    cout<<"energies1:  "<<den_en1.real()<<endl;
                    dcomp den_sps1=SP(vector_minus(ks1,ls),vector_minus(ks1,ls))*SP(vector_minus(ks1,q_in),vector_minus(ks1,q_in));
                    cout<<"sps1  :  "<<den_sps1.real()<<endl;
                    cout<<"num2:  "<<NLO_DT_DrellYan_numerator(ks2, ls, q_in)<<endl;
                    dcomp den_en2=2.*e_l*2.*(-e_l+q_in[0])*2.*(e_kl);
                    cout<<"energies2:  "<<den_en2.real()<<endl;
                    dcomp den_sps2=SP(ks2,ks2)*SP(vector_minus(ks2,q_in),vector_minus(ks2,q_in));
                    cout<<"sps2  :  "<<den_sps2.real()<<endl;
                    cout<<"num3:  "<<NLO_DT_DrellYan_numerator(ks3,ls,q_in)<<endl;
                    dcomp den_en3=2.*e_l*2.*(-e_l+q_in[0])*2.*(e_kq);
                    cout<<"energies3:  "<<den_en3.real()<<endl;
                    dcomp den_sps3=SP(ks3,ks3);
                    dcomp den_sps3p=SP(vector_minus(ks3,ls),vector_minus(ks3,ls));
                    cout<<"sps3  :  "<<den_sps3.real()<<endl;
                    cout<<e_l-e_kq-e_kl<<endl;
                    cout<<"sps3p  :  "<<den_sps3p.real()<<endl;
                    cout<<"res:  "<<rr.real()<<endl;
                    cout<<"-------"<<endl;
                };

                if(phase=0){
                    integrand_val=rr.real();
                }else{
                    integrand_val=rr.imag();
                };

                //MISSING DEFO JACOBIAN
                /**defo_jacques*/

                vector<dreal> out_obs_res;

                out_obs_res.push_back(integrand_val);
                out_obs_res.push_back(jac_flow*flow_h);
                out_obs_res.push_back(0);
                out_obs_res.push_back(0);
                out_obs_res.push_back(0);
                out_obs_res.push_back(0);

                return out_obs_res;

            };

            if(comp==3 || comp==4){

                vector<dreal> l_in=l;
                vector<dreal> k_in=k;
                vector<dreal> q_in={q0,0,0,0};

                if(comp==3){
                    l_in=l;
                    k_in=k;
                }else{
                    l_in=vector_swap(k);
                    k_in=vector_swap(l);
                    //q_in=vector_swap(q_in);
                };

                vector<vector<dreal>> ni={k_in,vector_swap(vector_minus(k_in,l_in)),vector_swap(l_in)};

                if(DEBUG==1){
                    cout<<"----"<<endl;
                    cout<<"presum: "<<norm(l_in)+norm(k_in)+norm(vector_minus(k_in,l_in))<<endl;
                };

                //vector<vector<dreal>> nf={};

                /*if(my_flow.t_val(ni,nf,q_in)==999){
                    vector<dcomp> res={0,0,0,0,0,0};
                    return res;
                };*/
                
                vector<dreal> ks=my_flow.perform_flow_jj(ni,q0,k_in);
                vector<dreal> ls=my_flow.perform_flow_jj(ni,q0,l_in);

                if(DEBUG==1){
                    cout<<ks[0]<<" "<<ks[1]<<" "<<ks[2]<<endl;
                    cout<<ls[0]<<" "<<ls[1]<<" "<<ls[2]<<endl;
                    cout<<ks[0]-ls[0]<<" "<<ks[1]-ls[1]<<" "<<ks[2]-ls[2]<<endl;
                };
                

                dreal jac_flow=my_flow.jacques_jj(ni,q0,2);
                dreal flow_h=my_flow.h(my_flow.t_val_jj(ni,q0));

                /*cout<<"t value:  "<<my_flow.t_val_jj(ni,q0)<<endl;
                cout<<"t jac value:  "<<my_flow.jacques_jj(ni,q0,2)<<endl;*/

                dreal e_k=norm(ks);
                dreal e_lq=norm(ls);
                dreal e_kl=norm(vector_minus(ls,ks));
                
                ks.insert(ks.begin(), e_k);
                ls.insert(ls.begin(), -e_lq+q_in[0]);


                vector<dreal> ks1={e_k,ks[0],ks[1],ks[2]};
                
                /*vector<vector<dreal>> constituents={ks,vector_swap(vector_minus(ks,ls)),vector_swap(vector_minus(ls,qs))};

                int spins[3]={1,0,-1};

                vector<dreal> obs_res=my_obs.b2b_sampling(constituents, 3, spins, a);*/

                if(DEBUG==1){
                    vector<dreal> ksp={0.3, 0.1, 0.5};
                    vector<dreal> lsp={0.2, 0.7, 0.1};
                    vector<dreal> qsp={0.1, 0.4, 0.5};

                    dreal ekp=norm(ksp);
                    dreal eqp=sqrt(pow(norm(qsp),2.)+pow(m,2));
                    dreal elp=norm(lsp);  
                    dreal ekqp=norm(vector_plus(ksp,qsp));
                    dreal elqp=norm(vector_minus(lsp,qsp)); //HERE  
                    dreal eklqp=norm(vector_plus(lsp,vector_plus(ksp,qsp)));

                    ksp.insert(ksp.begin(), ekp);
                    lsp.insert(lsp.begin(), -elqp+eqp);
                    qsp.insert(qsp.begin(), eqp);


                    cout<<"numerator: "<<NLO_DT_DrellYan_numerator(ksp,lsp,qsp)<<endl;
                    cout<<"energies: "<<2.*ekp*2.*(eqp-elqp-ekp)*2.*(elqp)*2*eqp<<endl;
                    cout<<"sps: "<<SP(vector_minus(ksp,qsp),vector_minus(ksp,qsp))*SP(lsp,lsp)<<endl;
                    cout<<NLO_DT_DrellYan_numerator(ksp,lsp,qsp)/
                    ((2.*ekp*2.*(eqp-elqp-ekp)*2.*(elp)*2*eqp)*
                    SP(vector_minus(ksp,qsp),vector_minus(ksp,qsp))*SP(lsp,lsp))<<endl; 

                };   

                dreal integrand_val;
                dcomp rr=NLO_DT_DrellYan_numerator(ks,ls,q_in)/
                    ((2.*e_k*2.*(q_in[0]-e_lq-e_k)*2.*(e_lq))*
                    SP(vector_minus(ks,q_in),vector_minus(ks,q_in))*SP(ls,ls)
                    );

                if(DEBUG==1){
                    cout<<"############REAL EMISSIONS########"<<endl;
                    cout<<"q0: "<<q0<<endl;
                    cout<<e_k+e_kl+e_lq<<endl;
                    cout<<"t: "<<my_flow.t_val_jj(ni,q0)<<endl;
                    cout<<"num:  "<<NLO_DT_DrellYan_numerator(ks,ls,q_in)<<endl;
                    dcomp den_en=2.*e_k*2.*(q_in[0]-e_lq-e_k)*2.*(e_lq);
                    cout<<"energies:  "<<den_en.real()<<endl;
                    dcomp den_sps=SP(vector_minus(ks,q_in),vector_minus(ks,q_in))*SP(ls,ls);
                    cout<<"sps  :  "<<den_sps.real()<<endl;
                    cout<<SP(vector_minus(ks,q_in),vector_minus(ks,q_in))<<endl;
                    cout<<SP(ls,ls)<<endl;
                    cout<<"res:  "<<rr.real()<<endl;
                    cout<<"-------"<<endl;
                };


                if(phase=0){
                    integrand_val=rr.real();
                }else{
                    integrand_val=rr.imag();
                };

                vector<dreal> out_obs_res;

                out_obs_res.push_back(integrand_val);
                out_obs_res.push_back(jac_flow*flow_h);
                out_obs_res.push_back(0);
                out_obs_res.push_back(0);
                out_obs_res.push_back(0);
                out_obs_res.push_back(0);

                return out_obs_res;

            };

        };


        observable_c NLO_se_DrellYan_integrand(vector<dreal> k, vector<dreal> l, vector<dreal> q, dreal a, int comp){

            if(comp==1){

                vector<dreal> l_in=l;
                vector<dreal> k_in=k;
                vector<dreal> q_in=q;

                dreal defo_sign;
                dreal jac_sign;

                if(comp==1){
                    l_in=l;
                    k_in=k;
                    q_in=q;
                    defo_sign=1;
                    jac_sign=1;
                }else{
                    l_in=vector_swap(k);
                    k_in=vector_swap(l);
                    q_in=vector_swap(q);
                    defo_sign=1;
                    jac_sign=-1;
                };

                vector<vector<dreal>> ni={k_in,vector_swap(vector_minus(k_in,q_in))};
                vector<vector<dreal>> nf={};

                if(my_flow.t_val(ni,nf,q_in)==999){
                    observable_c res;
                    vector<dreal> j1={0,0,0,0};
                    vector<dreal> j2={0,0,0,0};
                    vector<dreal> pg={0,0,0,0};
                    res.eval=0;
                    res.jac=0;
                    res.j1=j1;
                    res.j2=j2;
                    res.spin1=0;
                    res.spin2=0;
                    res.pg=pg;
                    return res;
                };
                
                vector<dreal> ks=my_flow.perform_flow(ni,nf,q_in,k_in);
                vector<dreal> ls_real=my_flow.perform_flow(ni,nf,q_in,l_in);
                vector<dreal> qs=my_flow.perform_flow(ni,nf,q_in,q_in);

                dreal jac_flow=my_flow.jacques(ni,nf,q_in,3);
                dreal flow_h=my_flow.h(my_flow.t_val(ni,nf,q_in));

                dcomp i(0,1.);
                vector<dreal> kappa={my_defo.get_deformation_dt(ks,ls_real,qs,1,1).rpart(),my_defo.get_deformation_dt(ks,ls_real,qs,1,2).rpart(),my_defo.get_deformation_dt(ks,ls_real,qs,1,3).rpart()};
                vector<dcomp> ls={ls_real[0]+defo_sign*i*kappa[0],ls_real[1]+defo_sign*i*kappa[1],ls_real[2]+defo_sign*i*kappa[2]};
                dcomp defo_jacques=my_defo.jacques_deformation_dt(ls_real,ks,qs,1);

                dreal e_k=norm(ks);
                dcomp e_kq=norm(vector_minus(ks,qs));
                dcomp e_l=norm(ls);
                dcomp e_kl=norm(vector_minus(ls,ks));
                dreal e_q=sqrt(pow(norm(qs),2.)+pow(m,2));

                dcomp e_kUV=Sqrt(pow(norm(ks),2.)+pow(MUV,2));
                
                
                ks.insert(ks.begin(), e_k);
                qs.insert(qs.begin(), e_q);

                vector<dcomp> ls1={e_l,ls[0],ls[1],ls[2]};

                dcomp se_virtual_1=-NLO_SE_DrellYan_numerator(ks,ls1,qs)*2.*(ks[0]-ls1[0])/
                    (2.*e_l*2.*e_k*2.*e_k*2.*(-e_k+e_q)*2.*e_q*SP(vector_minus(ls1,ks),vector_minus(ls1,ks))*SP(vector_minus(ls1,ks),vector_minus(ls1,ks)));

                vector<dcomp> ls2={e_k+e_kl,ls[0],ls[1],ls[2]};

                dcomp se_virtual_2=-1/(2.*e_k*2.*e_k*2.*(-e_k+e_q)*2.*e_q)*(NLO_SE_DrellYan_numerator(ks, ls2, qs)*2.*ls2[0]/
                    (2.*e_kl*SP(ls2,ls2)*SP(ls2,ls2))-NLO_SE_DrellYan_numerator_d(ks, ls2, qs)/(2.*e_kl*SP(ls2,ls2)));

                dcomp e_lUV=Sqrt(pow(norm(ls),2.)+pow(MUV,2.));
                vector<dcomp> lsuv={e_lUV,ls[0],ls[1],ls[2]};

                dcomp se_uv=1/(8.*e_k*e_k*(-e_k+e_q)*2*e_q)*(
                            -2.*(-e_k+e_lUV)*(48.*e_k*SP(ks,qs)/Power(e_lUV,4.)-48*SP(ks,qs)*SP(ks,lsuv)/Power(e_lUV,5.))  
                            -4.*(-16.*e_k*SP(ks,qs)/Power(e_lUV,3.)+24.*SP(ks,qs)*SP(ks,lsuv)/Power(e_lUV,4.)));

                
                //cout<<0.5*pow(e_l,3.)*se_uv<<" "<<pow(e_l,3.)*(se_virtual_1+se_virtual_2)<<endl;

                vector<vector<dreal>> constituents={ks,vector_swap(vector_minus(ks,qs))};

                int spins[2]={1,-1};

                vector<dreal> pg={0,0,0,0};

                observable obs_res_n=my_obs.b2b_sampling_final(constituents, 2, pg, spins, a);

                observable_c obs_res;
                obs_res.eval=obs_res_n.eval;
                obs_res.jac=obs_res_n.jac;
                obs_res.j1=obs_res_n.j1;
                obs_res.j2=obs_res_n.j2;
                obs_res.spin1=obs_res_n.spin1;
                obs_res.spin2=obs_res_n.spin2;
                obs_res.pg=obs_res_n.pg;
                
                if(DEBUG==1){
                    vector<dreal> ksp={0.3, 0.1, 0.5};
                    vector<dreal> lsp={0.2, 0.7, 0.1};
                    vector<dreal> qsp={0.4, 0.4, 0.1};

                    dreal ekp=norm(ksp);
                    dreal ekqp=norm(vector_minus(ksp,qsp));
                    dreal eqp=sqrt(pow(norm(qsp),2.)+pow(m,2));
                    dreal elp=norm(lsp);  
                    dreal eklp=norm(vector_minus(ksp,lsp));  
                    dreal eklqp=norm(vector_plus(lsp,vector_plus(ksp,qsp)));
                    dreal e_kUVp=sqrt(pow(norm(ksp),2.)+pow(MUV,2));
                    
                    ksp.insert(ksp.begin(), ekp);
                    qsp.insert(qsp.begin(), eqp);
                    

                    vector<dreal> ls1p={elp,lsp[0],lsp[1],lsp[2]};

                    dreal se_virtual_1p=-NLO_SE_DrellYan_numerator(ksp,ls1p,qsp)*2.*(ksp[0]-ls1p[0])/
                        (2.*elp*2.*ekp*2.*ekp*2.*(-ekp+eqp)*2.*eqp*SP(vector_minus(ls1p,ksp),vector_minus(ls1p,ksp))*SP(vector_minus(ls1p,ksp),vector_minus(ls1p,ksp)));

                    vector<dreal> ls2p={ekp+eklp,lsp[0],lsp[1],lsp[2]};

                    dreal se_virtual_2p=-1/(2.*ekp*2.*ekp*2.*(-ekp+eqp)*2.*eqp)*(NLO_SE_DrellYan_numerator(ksp, ls2p, qsp)*2.*ls2p[0]/
                        (2.*eklp*SP(ls2p,ls2p)*SP(ls2p,ls2p))-NLO_SE_DrellYan_numerator_d(ksp, ls2p, qsp)/(2.*eklp*SP(ls2p,ls2p)));

                    dreal e_lUVp=Sqrt(pow(norm(lsp),2.)+pow(MUV,2.));
                    vector<dreal> lsuvp={e_lUVp,lsp[0],lsp[1],lsp[2]};
                    

                    cout<<"ltd1: "<<se_virtual_1p<<endl;
                    cout<<NLO_SE_DrellYan_numerator(ksp,ls1p,qsp)<<endl;
                    cout<<2.*elp*2.*ekp*2.*ekp*2.*(-ekp+eqp)*2.*eqp<<endl;
                    cout<<SP(vector_minus(ls1p,ksp),vector_minus(ls1p,ksp))*SP(vector_minus(ls1p,ksp),vector_minus(ls1p,ksp))<<endl;
                    cout<<"ltd2: "<<se_virtual_2p<<endl;
                    cout<<1/(2.*ekp*2.*ekp*2.*(-ekp+eqp)*2.*eqp)*NLO_SE_DrellYan_numerator_d(ksp, ls2p, qsp)/(2.*eklp*SP(ls2p,ls2p))<<endl;
                    cout<<-1/(2.*ekp*2.*ekp*2.*(-ekp+eqp)*2.*eqp)*NLO_SE_DrellYan_numerator(ksp, ls2p, qsp)*2.*ls2p[0]/(2.*eklp*SP(ls2p,ls2p)*SP(ls2p,ls2p))<<endl;
                    cout<<NLO_SE_DrellYan_numerator(ksp, ls2p, qsp)<<endl;
                    cout<<"res: "<<se_virtual_1p+se_virtual_2p<<endl;
                    cout<<e_lUVp<<endl;
                    cout<<1/(8.*ekp*ekp*(-ekp+eqp)*2*eqp)<<endl;
                    cout<<(3./(32.*Power(e_lUVp,4.))-3.*(-ekp+e_lUVp)/(32.*Power(e_lUVp,5.)))<<endl;
                    cout<<1/(8.*ekp*ekp*(-ekp+eqp)*2*eqp)*(
                            -2.*(-ekp+e_lUVp)*(48.*ekp*SP(ksp,qsp)/Power(e_lUVp,4.)-48*SP(ksp,qsp)*SP(ksp,lsuvp)/Power(e_lUVp,5.))  
                            -4.*(-16.*ekp*SP(ksp,qsp)/Power(e_lUVp,3.)+24.*SP(ksp,qsp)*SP(ksp,lsuvp)/Power(e_lUVp,4.)))<<endl;
                    cout<<-2.*(-ekp+e_lUVp)<<endl;
                    cout<<48.*ekp*SP(ksp,qsp)/Power(e_lUVp,4.)<<endl;
                    cout<<48*SP(ksp,qsp)*SP(ksp,lsuv)/Power(e_lUVp,5.)<<endl;
                    cout<<-2.*(-ekp+e_lUVp)*(48.*ekp*SP(ksp,qsp)/Power(e_lUVp,4.)-48*SP(ksp,qsp)*SP(ksp,lsuv)/Power(e_lUVp,5.))<<endl;
                    cout<<-4.*(-16.*ekp*SP(ksp,qsp)/Power(e_lUVp,3.)+24.*SP(ksp,qsp)*SP(ksp,lsuv)/Power(e_lUVp,4.))<<endl;

                };

                dreal colour_average=group_average(obs_res.spin1,obs_res.spin2);

                dcomp integrand_val=constant*colour_average*(se_virtual_1+se_virtual_2+0.5*se_uv);

                obs_res.jac=obs_res.jac*jac_flow*flow_h*defo_jacques;

                obs_res.eval=integrand_val;

                return obs_res;


            };


            if(comp==2){

                vector<dreal> l_in=l;
                vector<dreal> k_in=k;
                vector<dreal> q_in=q;

                vector<vector<dreal>> ni={l_in,vector_minus(k_in,l_in),vector_swap(vector_minus(k_in,q_in))};
                vector<vector<dreal>> nf={};

                if(my_flow.t_val(ni,nf,q_in)==999){
                    observable_c res;
                    vector<dreal> j1={0,0,0,0};
                    vector<dreal> j2={0,0,0,0};
                    vector<dreal> pg={0,0,0,0};
                    res.eval=0;
                    res.jac=0;
                    res.j1=j1;
                    res.j2=j2;
                    res.spin1=0;
                    res.spin2=0;
                    res.pg=pg;
                    return res;
                };
                
                vector<dreal> ks=my_flow.perform_flow(ni,nf,q_in,k_in);
                vector<dreal> ls=my_flow.perform_flow(ni,nf,q_in,l_in);
                vector<dreal> qs=my_flow.perform_flow(ni,nf,q_in,q_in);

                dreal jac_flow=my_flow.jacques(ni,nf,q_in,3);
                dreal flow_h=my_flow.h(my_flow.t_val(ni,nf,q_in));

                dreal e_k=norm(ks);
                dreal e_kq=norm(vector_minus(ks,qs));
                dreal e_l=norm(ls);
                dreal e_kl=norm(vector_minus(ls,ks));
                dreal e_q=sqrt(pow(norm(qs),2.)+pow(m,2));

                dcomp e_kUV=Sqrt(pow(norm(ks),2.)+pow(MUV,2));
                
                
                ks.insert(ks.begin(), e_l+e_kl);
                ls.insert(ls.begin(), e_l);
                qs.insert(qs.begin(), e_q);

                dcomp integrand_val=NLO_SE_DrellYan_numerator(ks,ls,qs)/(2.*e_l*2.*e_kl*2.*(-e_kl-e_l+e_q)*2*e_q*SP(ks,ks)*SP(ks,ks)); 
                
                vector<vector<dreal>> constituents={ls,vector_minus(ks,ls),vector_swap(vector_minus(ks,qs))};

                int spins[3]={1,0,-1};

                vector<dreal> pg={0,0,0,0};

                observable obs_res_n=my_obs.b2b_sampling_final(constituents, 3, pg, spins, a);


                observable_c obs_res;
                obs_res.eval=obs_res_n.eval;
                obs_res.jac=obs_res_n.jac;
                obs_res.j1=obs_res_n.j1;
                obs_res.j2=obs_res_n.j2;
                obs_res.spin1=obs_res_n.spin1;
                obs_res.spin2=obs_res_n.spin2;
                obs_res.pg=obs_res_n.pg;
                
                if(DEBUG==1){
                    vector<dreal> ksp={0.3, 0.1, 0.5};
                    vector<dreal> lsp={0.2, 0.7, 0.1};
                    vector<dreal> qsp={0.4, 0.4, 0.1};

                    dreal ekp=norm(ksp);
                    dreal ekqp=norm(vector_minus(ksp,qsp));
                    dreal eqp=sqrt(pow(norm(qsp),2.)+pow(m,2));
                    dreal elp=norm(lsp);  
                    dreal eklp=norm(vector_minus(ksp,lsp));  
                    dreal eklqp=norm(vector_plus(lsp,vector_plus(ksp,qsp)));
                    dreal e_kUVp=sqrt(pow(norm(ksp),2.)+pow(MUV,2));

                    ksp.insert(ksp.begin(), elp+eklp);
                    lsp.insert(lsp.begin(), elp);
                    qsp.insert(qsp.begin(), eqp);

                    dreal res=NLO_SE_DrellYan_numerator(ksp,lsp,qsp)/(2.*elp*2.*eklp*2.*(-eklp-elp+eqp)*SP(ksp,ksp)*SP(ksp,ksp));

                    cout<<res<<endl;

                };

                obs_res.jac=obs_res.jac*jac_flow*flow_h;

                dreal colour_average=group_average(obs_res.spin1,obs_res.spin2);

                obs_res.eval=constant*colour_average*integrand_val;

                return obs_res;


            };

            observable_c obs_res_f;
            obs_res_f.eval=0;
            obs_res_f.jac=0;
            obs_res_f.j1={0,0,0,0};
            obs_res_f.j2={0,0,0,0};
            obs_res_f.spin1=0;
            obs_res_f.spin2=0;
            obs_res_f.pg={0,0,0,0};

            return obs_res_f;

        };


        vector<dreal> NLO_se_jj_integrand(vector<dreal> k, vector<dreal> l, dreal q0, int comp, int phase){

            if(comp==1){

                vector<dreal> l_in=l;
                vector<dreal> k_in=k;
                vector<dreal> q_in={q0,0,0,0};

                dreal defo_sign;

                vector<vector<dreal>> ni={k_in,vector_swap(k_in)};
                
                
                vector<dreal> ks=my_flow.perform_flow_jj(ni,q0,k_in);
                vector<dreal> ls=my_flow.perform_flow_jj(ni,q0,l_in);

                dreal jac_flow=my_flow.jacques_jj(ni,q0,2);
                dreal flow_h=my_flow.h(my_flow.t_val_jj(ni,q0));

                /*dcomp i(0,1.);
                vector<dreal> kappa=my_defo.get_deformation(ls_real,ks,1);
                vector<dcomp> ls={ls_real[0]+defo_sign*i*kappa[0],ls_real[1]+defo_sign*i*kappa[1],ls_real[2]+defo_sign*i*kappa[2]};
                dcomp defo_jacques=my_defo.jacques_deformation(ls_real,ks,1,defo_sign);*/

                dreal e_k=norm(ks);
                dcomp e_kq=norm(ks);
                dcomp e_l=norm(ls);
                dcomp e_kl=norm(vector_minus(ls,ks));
                

                dcomp e_kUV=Sqrt(pow(norm(ks),2.)+pow(MUV,2));
                
                
                ks.insert(ks.begin(), e_k);
                

                vector<dcomp> ls1={e_l,ls[0],ls[1],ls[2]};

                dcomp se_virtual_1=-NLO_SE_DrellYan_numerator(ks,ls1,q_in)*2.*(ks[0]-ls1[0])/
                    (2.*e_l*2.*e_k*2.*e_k*2.*(-e_k+q_in[0])*SP(vector_minus(ls1,ks),vector_minus(ls1,ks))*SP(vector_minus(ls1,ks),vector_minus(ls1,ks)));

                vector<dcomp> ls2={e_k+e_kl,ls[0],ls[1],ls[2]};

                dcomp se_virtual_2=-1/(2.*e_k*2.*e_k*2.*(-e_k+q_in[0]))*(NLO_SE_DrellYan_numerator(ks, ls2, q_in)*2.*ls2[0]/
                    (2.*e_kl*SP(ls2,ls2)*SP(ls2,ls2))-NLO_SE_DrellYan_numerator_d(ks, ls2, q_in)/(2.*e_kl*SP(ls2,ls2)));

                dcomp e_lUV=Sqrt(pow(norm(ls),2.)+pow(MUV,2.));
                vector<dcomp> lsuv={e_lUV,ls[0],ls[1],ls[2]};
                //vector<dcomp> lsuv={0,ls[0],ls[1],ls[2]};

                dcomp se_uv=(1/2.)*1/(8.*e_k*e_k*(-e_k+q_in[0]))*(
                            -2.*(-e_k+e_lUV)*(48.*e_k*SP(ks,q_in)/Power(e_lUV,4.)-48*SP(ks,q_in)*SP(ks,lsuv)/Power(e_lUV,5.))  
                            -4.*(-16.*e_k*SP(ks,q_in)/Power(e_lUV,3.)+24.*SP(ks,q_in)*SP(ks,lsuv)/Power(e_lUV,4.)));
                
                if(DEBUG==1){
                    cout<<"virtual:  "<<se_virtual_1.real()+se_virtual_2.real()<<endl;
                    cout<<"uv virtual:  "<<se_uv.real()<<endl;
                    cout<<"uv virtual1:  "<<1/(8.*e_k*e_k*(-e_k+q_in[0]))*(
                                -2.*(-e_k+e_lUV)*(48.*e_k*SP(ks,q_in)/Power(e_lUV,4.)))<<endl;
                    cout<<"uv virtual2:  "<<1/(8.*e_k*e_k*(-e_k+q_in[0]))*(
                                -2.*(-e_k+e_lUV)*(-48*SP(ks,q_in)*SP(ks,lsuv)/Power(e_lUV,5.)))<<endl;
                    cout<<"uv virtual3:  "<<1/(8.*e_k*e_k*(-e_k+q_in[0]))*(
                                -4.*(-16.*e_k*SP(ks,q_in)/Power(e_lUV,3.)))<<endl;
                    cout<<"uv virtual4:  "<<1/(8.*e_k*e_k*(-e_k+q_in[0]))*(
                                -4.*(24.*SP(ks,q_in)*SP(ks,lsuv)/Power(e_lUV,4.)))<<endl;
                    cout<<"jacobian:  "<<jac_flow*flow_h<<endl;
                };
                


                dcomp integrand_val=se_virtual_1+se_virtual_2+se_uv;

                dreal rr;

                if(phase==0){
                    rr=integrand_val.real();
                }else{
                    rr=integrand_val.imag();
                };


                vector<dreal> out_obs_res;

                out_obs_res.push_back(rr);
                out_obs_res.push_back(jac_flow*flow_h);
                out_obs_res.push_back(0);
                out_obs_res.push_back(0);
                out_obs_res.push_back(0);
                out_obs_res.push_back(0);

                return out_obs_res;


            };


            if(comp==2){

                vector<dreal> l_in=l;
                vector<dreal> k_in=k;
                vector<dreal> q_in={q0,0,0,0};

                vector<vector<dreal>> ni={l_in,vector_minus(k_in,l_in),vector_swap(k_in)};
                
                vector<dreal> ks=my_flow.perform_flow_jj(ni,q0,k_in);
                vector<dreal> ls=my_flow.perform_flow_jj(ni,q0,l_in);
   

                dreal jac_flow=my_flow.jacques_jj(ni,q0,2);
                dreal flow_h=my_flow.h(my_flow.t_val_jj(ni,q0));

                dreal e_k=norm(ks);
                dreal e_kq=norm(ks);
                dreal e_l=norm(ls);
                dreal e_kl=norm(vector_minus(ls,ks));

                dcomp e_kUV=Sqrt(pow(norm(ks),2.)+pow(MUV,2));
                
                
                ks.insert(ks.begin(), e_l+e_kl);
                ls.insert(ls.begin(), e_l);
  

                dcomp integrand_val=NLO_SE_DrellYan_numerator(ks,ls,q_in)/(2.*e_l*2.*e_kl*2.*(-e_kl-e_l+q_in[0])*SP(ks,ks)*SP(ks,ks)); 


                dreal rr;

                if(phase==0){
                    rr=integrand_val.real();
                }else{
                    rr=integrand_val.imag();
                };

                vector<dreal> out_obs_res;

                out_obs_res.push_back(rr);
                out_obs_res.push_back(jac_flow*flow_h);
                out_obs_res.push_back(0);
                out_obs_res.push_back(0);
                out_obs_res.push_back(0);
                out_obs_res.push_back(0);

                return out_obs_res;


            };



        };

        /*
        vector<dreal> NLO_se_DrellYan_integrand(vector<dreal> k, vector<dreal> l, vector<dreal> q, dreal a, int comp){

            if(comp==1){

                vector<vector<dreal>> ni={k,vector_swap(vector_minus(k,q))};
                vector<vector<dreal>> nf={};

                if(my_flow.t_val(ni,nf,q)==999){
                    vector<dreal> res={0,0,0,0,0,0};
                    return res;
                };
                
                vector<dreal> ks=my_flow.perform_flow(ni,nf,q,k);
                vector<dreal> ls=my_flow.perform_flow(ni,nf,q,l);
                vector<dreal> qs=my_flow.perform_flow(ni,nf,q,q);

                dreal jac_flow=my_flow.jacques(ni,nf,q,3);
                dreal flow_h=my_flow.h(my_flow.t_val(ni,nf,q));

                dreal e_k=norm(ks);
                dreal e_l=norm(ls);
                dreal e_kl=norm(vector_plus(ls,ks));
                dreal e_q=sqrt(pow(norm(qs),2.)+pow(m,2));
                
                ks.insert(ks.begin(), e_k);
                qs.insert(ls.begin(), e_q);

                vector<dreal> ls1={e_k,ks[0],ks[1],ks[2]};


                
                vector<vector<dreal>> constituents={ks,vector_swap(vector_minus(ks,qs))};

                int spins[2]={1,-1};

                vector<dreal> obs_res=my_obs.b2b_sampling(constituents, 2, spins, a);

                if(DEBUG==1){


                };   


                    
                dreal integrand_val=NLO_DT_DrellYan_numerator(ks,ls,qs)/
                    ((2.*e_k*2.*(e_q-e_l-e_k)*2.*(e_l)*2*e_q)*
                    SP(vector_minus(ks,qs),vector_minus(ks,qs))*SP(ls,ls)
                    );

                obs_res[0]=obs_res[0]*jac_flow*flow_h;

                obs_res.insert(obs_res.begin(), 
                    integrand_val
                    );

                return obs_res;

            };

        };*/



        /*vector<dreal> NLO_dt_DrellYan_integrand*/

        observable NLO_st_channel_DrellYan_integrand(vector<dreal> k, vector<dreal> l, vector<dreal> q, dreal a, int comp){
            
            if(comp==1){
                
                vector<vector<dreal>> ni={vector_swap(k),vector_plus(vector_plus(k,l),q)};
                vector<vector<dreal>> nf={l};

                if(my_flow.t_val(ni,nf,q)==999){
                    observable res;
                    vector<dreal> j1={0,0,0,0};
                    vector<dreal> j2={0,0,0,0};
                    vector<dreal> pg={0,0,0,0};
                    res.eval=0;
                    res.jac=0;
                    res.j1=j1;
                    res.j2=j2;
                    res.spin1=0;
                    res.spin2=0;
                    res.pg=pg;
                    return res;
                };

                vector<dreal> ks=my_flow.perform_flow(ni,nf,q,k);
                vector<dreal> ls=my_flow.perform_flow(ni,nf,q,l);
                vector<dreal> qs=my_flow.perform_flow(ni,nf,q,q);

                dreal jac_flow=my_flow.jacques(ni,nf,q,3);
                dreal flow_h=my_flow.h(my_flow.t_val(ni,nf,q));

                if(DEBUG==1){
                    cout<<"k:  "<<k[0]<<" "<<k[1]<<" "<<k[2]<<endl;
                    cout<<"l:  "<<l[0]<<" "<<l[1]<<" "<<l[2]<<endl;
                    cout<<"q:  "<<q[0]<<" "<<q[1]<<" "<<q[2]<<endl;
                    cout<<"scaled_momenta"<<endl;
                    cout<<"k:  "<<ks[0]<<" "<<ks[1]<<" "<<ks[2]<<endl;
                    cout<<"l:  "<<ls[0]<<" "<<ls[1]<<" "<<ls[2]<<endl;
                    cout<<"q:  "<<qs[0]<<" "<<qs[1]<<" "<<qs[2]<<endl;
                    cout<<"t_value"<<endl;
                    cout<<my_flow.t_val(ni,nf,q)<<endl;
                    cout<<"scaling jacobian"<<endl;
                    cout<<my_flow.jacques(ni,nf,q,3)<<endl;
                };

                dreal ek=norm(ks);
                dreal eq=sqrt(pow(norm(qs),2.)+pow(m,2));
                dreal el=norm(ls);  

                ks.insert(ks.begin(), -ek);
                ls.insert(ls.begin(), el);
                qs.insert(qs.begin(), eq);

                vector<vector<dreal>> constituents={vector_swap(ks),vector_plus(vector_plus(ks,ls),qs)};

                int spins[2]={1,0};

                observable obs_res=my_obs.b2b_sampling_final(constituents, 2, ls, spins, a);;

                if(DEBUG==1){
                    vector<dreal> ksp={0.3, 0.1, 0.5};
                    vector<dreal> lsp={0.2, 0.7, 0.1};
                    vector<dreal> qsp={0.1, 0.4, 0.5};

                    dreal ekp=norm(ksp);
                    dreal eqp=sqrt(pow(norm(qsp),2.)+pow(1,2));
                    dreal elp=norm(lsp);  
                    dreal ekqp=norm(vector_plus(ksp,qsp));
                    dreal elqp=norm(vector_plus(lsp,qsp));  
                    dreal eklqp=norm(vector_plus(lsp,vector_plus(ksp,qsp)));

                    ksp.insert(ksp.begin(), -ekp);
                    lsp.insert(lsp.begin(), elp);
                    qsp.insert(qsp.begin(), eqp);

                    cout<<"numerator: "<<NLO_DT_DrellYan_numerator(ksp,vector_swap(vector_plus(lsp,qsp)),vector_swap(qsp))<<endl;

                    cout<<NLO_DT_DrellYan_numerator(ksp,vector_swap(vector_plus(lsp,qsp)),vector_swap(qsp))/(8*ekp*elp*(-ekp+elp+eqp)*eqp*SP(vector_plus(ksp,qsp),vector_plus(ksp,qsp))*SP(vector_plus(lsp,qsp),vector_plus(lsp,qsp)))<<endl;
                };
                
                dreal integrand_val=constant*NLO_DT_DrellYan_numerator(ks,vector_swap(vector_plus(ls,qs)),vector_swap(qs))/
                    (16*abs(ek*el*(-ek+el+eq)*eq)*SP(vector_plus(ks,qs),vector_plus(ks,qs))*SP(vector_plus(ls,qs),vector_plus(ls,qs)));
                
                //cout<<"ks+qs: "<<vector_plus(ks,qs)[0]<<" "<<vector_plus(ks,qs)[1]<<" "<<vector_plus(ks,qs)[2]<<" "<<vector_plus(ks,qs)[3]<<endl;
                //cout<<"sp: "<<SP(vector_plus(ks,qs),vector_plus(ks,qs))<<endl;
                //cout<<"integrand val: "<<integrand_val<<endl;

                if(DEBUG==1){
                    vector<dreal> ksp={-121.702, 0., 0., -121.702};
                    vector<dreal> lsp={104.621, -22.4313, -15.3108, -101.035};
                    vector<dreal> qsp={138.784, 22.4313, 15.3108, 101.035};

                    cout<<"comparison_MG: "<<constant*NLO_DT_DrellYan_numerator(ksp,vector_swap(vector_plus(lsp,qsp)),vector_swap(qsp))/(SP(vector_plus(ksp,qsp),vector_plus(ksp,qsp))*SP(vector_plus(lsp,qsp),vector_plus(lsp,qsp)))<<endl;
                };

                if(DEBUG==1){
                    vector<dreal> ksp={-70.962, 0., 0., -70.962};
                    vector<dreal> lsp={41.6672, 37.4207, -18.0886, -2.94203};
                    vector<dreal> qsp={100.25671940723308, -37.4207, 18.0886, 2.94203};

                    cout<<"comparison_MG_2: "<<constant*NLO_DT_DrellYan_numerator(ksp,vector_swap(vector_plus(lsp,qsp)),vector_swap(qsp))/(SP(vector_plus(ksp,qsp),vector_plus(ksp,qsp))*SP(vector_plus(lsp,qsp),vector_plus(lsp,qsp)))<<endl;
                };

                if(DEBUG==1){
                    cout<<"############### cut 1 info #############"<<endl;
                    cout<<"original kinematics"<<endl;
                    cout<<"k:  "<<k[0]<<" "<<k[1]<<" "<<k[2]<<endl;
                    cout<<"l:  "<<l[0]<<" "<<l[1]<<" "<<l[2]<<endl;
                    cout<<"q:  "<<q[0]<<" "<<q[1]<<" "<<q[2]<<endl;
                    cout<<"kinematics"<<endl;
                    cout<<"k:  "<<ks[0]<<" "<<ks[1]<<" "<<ks[2]<<" "<<ks[3]<<endl;
                    cout<<"l:  "<<ls[0]<<" "<<ls[1]<<" "<<ls[2]<<" "<<ls[3]<<endl;
                    cout<<"q:  "<<qs[0]<<" "<<qs[1]<<" "<<qs[2]<<" "<<qs[3]<<endl;
                    cout<<"important momenta"<<endl;
                    cout<<"k:  "<<ks[0]+qs[0]<<" "<<ks[1]+qs[1]<<" "<<ks[2]+qs[2]<<" "<<ks[3]+qs[3]<<endl;
                    cout<<"l:  "<<ls[0]<<" "<<ls[1]<<" "<<ls[2]<<" "<<ls[3]<<endl;
                    cout<<"t value: "<<my_flow.t_val(ni,nf,q)<<endl;
                    cout<<"flow jacobian: "<<jac_flow<<endl;
                    cout<<"h function: "<<flow_h<<endl;
                    cout<<"back 2 back jacobian: "<<obs_res.jac<<endl;
                    cout<<"numerator"<<NLO_DT_DrellYan_numerator(ks,vector_swap(vector_plus(ls,qs)),vector_swap(qs))<<endl;
                    cout<<constant<<endl;
                    cout<<"integrand value: "<<integrand_val/constant<<endl;
                    cout<<"integrand value: "<<integrand_val<<endl;
                    cout<<"product of energies: "<<abs(ek*el*(-ek+el+eq)*eq)<<endl;
                    cout<<"product of propagators: "<<SP(vector_plus(ks,qs),vector_plus(ks,qs))*SP(vector_plus(ls,qs),vector_plus(ls,qs))<<endl;
                    cout<<"final result: "<<jac_flow*flow_h*obs_res.jac*integrand_val<<endl;
                };

                if(DEBUG==2){
                    cout<<"out clustering"<<endl;
                    cout<<obs_res.eval<<endl;
                    cout<<obs_res.jac<<endl;
                };

                obs_res.jac=obs_res.jac*jac_flow*flow_h;

                dreal colour_average=group_average(obs_res.spin1,obs_res.spin2);

                obs_res.eval=colour_average*integrand_val;

                if(DEBUG==2){
                    cout<<"final"<<endl;
                    cout<<colour_average<<endl;
                    cout<<obs_res.eval<<endl;
                    cout<<obs_res.jac<<endl;
                    cout<<obs_res.jac<<endl;
                    cout<<"------"<<endl;
                };

                

                return obs_res;

            };

            if(comp==2){
                vector<vector<dreal>> ni={vector_swap(k),vector_plus(k,q)};
                vector<vector<dreal>> nf={};

                if(my_flow.t_val(ni,nf,q)==999){
                    observable res;
                    vector<dreal> j1={0,0,0,0};
                    vector<dreal> j2={0,0,0,0};
                    vector<dreal> pg={0,0,0,0};
                    res.eval=0;
                    res.jac=0;
                    res.j1=j1;
                    res.j2=j2;
                    res.spin1=0;
                    res.spin2=0;
                    res.pg=pg;
                    return res;
                };

                vector<dreal> ks=my_flow.perform_flow(ni,nf,q,k);
                vector<dreal> ls=my_flow.perform_flow(ni,nf,q,l);
                vector<dreal> qs=my_flow.perform_flow(ni,nf,q,q);

                dreal jac_flow=my_flow.jacques(ni,nf,q,3);
                dreal flow_h=my_flow.h(my_flow.t_val(ni,nf,q));

                dreal ek=norm(ks);
                dreal eq=sqrt(pow(norm(qs),2.)+pow(m,2));
                dreal el=norm(ls);  

                ks.insert(ks.begin(), -ek);
                ls.insert(ls.begin(), el);
                qs.insert(qs.begin(), eq);

                vector<vector<dreal>> constituents={vector_swap(ks),ls,vector_plus(ks,qs)};

                int spins[3]={1,1,-1};

                observable obs_res=my_obs.b2b_sampling_final(constituents, 3, ls, spins, a);;

                if(DEBUG==1){
                    vector<dreal> ksp={0.3, 0.1, 0.5};
                    vector<dreal> lsp={0.2, 0.7, 0.1};
                    vector<dreal> qsp={0.1, 0.4, 0.5};

                    dreal ekp=norm(ksp);
                    dreal eqp=sqrt(pow(norm(qsp),2.)+pow(m,2));
                    dreal elp=norm(lsp);  
                    dreal ekqp=norm(vector_plus(ksp,qsp));
                    dreal elqp=norm(vector_plus(lsp,qsp));  
                    dreal eklqp=norm(vector_plus(lsp,vector_plus(ksp,qsp)));

                    ksp.insert(ksp.begin(), -ekp);
                    lsp.insert(lsp.begin(), elp);
                    qsp.insert(qsp.begin(), eqp);

                    cout<<NLO_DT_DrellYan_numerator(ksp,vector_swap(vector_plus(lsp,qsp)),vector_swap(qsp))/(16*ekp*(-ekp+eqp)*elp*eqp*SP(vector_plus(qsp,vector_plus(ksp,lsp)),vector_plus(qsp,vector_plus(ksp,lsp)))*SP(vector_plus(qsp,lsp),vector_plus(qsp,lsp)))<<endl;
                };
                //cout<<"hallo"<<endl;


                dreal integrand_val=constant*NLO_DT_DrellYan_numerator(ks,vector_swap(vector_plus(ls,qs)),vector_swap(qs))/(16*abs(ek*(-ek+eq)*el*eq)*SP(vector_plus(qs,vector_plus(ks,ls)),vector_plus(qs,vector_plus(ks,ls)))*SP(vector_plus(qs,ls),vector_plus(qs,ls)));

                if(DEBUG==1){
                    cout<<"############### cut 2 info #############"<<endl;
                    cout<<"original kinematics"<<endl;
                    cout<<"k:  "<<k[0]<<" "<<k[1]<<" "<<k[2]<<endl;
                    cout<<"l:  "<<l[0]<<" "<<l[1]<<" "<<l[2]<<endl;
                    cout<<"q:  "<<q[0]<<" "<<q[1]<<" "<<q[2]<<endl;
                    cout<<"rescaled kinematics"<<endl;
                    cout<<"k:  "<<ks[0]<<" "<<ks[1]<<" "<<ks[2]<<" "<<ks[3]<<endl;
                    cout<<"l:  "<<ls[0]<<" "<<ls[1]<<" "<<ls[2]<<" "<<ls[3]<<endl;
                    cout<<"q:  "<<qs[0]<<" "<<qs[1]<<" "<<qs[2]<<" "<<qs[3]<<endl;
                    cout<<"important momenta"<<endl;
                    cout<<"k+q:  "<<ks[0]+qs[0]<<" "<<ks[1]+qs[1]<<" "<<ks[2]+qs[2]<<" "<<ks[3]+qs[3]<<endl;
                    cout<<"l:  "<<ls[0]<<" "<<ls[1]<<" "<<ls[2]<<" "<<ls[3]<<endl;
                    cout<<"t value: "<<my_flow.t_val(ni,nf,q)<<endl;
                    cout<<"flow jacobian: "<<jac_flow<<endl;
                    cout<<"h function: "<<flow_h<<endl;
                    cout<<"back 2 back jacobian: "<<obs_res.jac<<endl;
                    cout<<"integrand value: "<<integrand_val<<endl;
                    cout<<"product of energies: "<<abs(ek*(-ek+eq)*el*eq)<<endl;
                    cout<<"product of propagators: "<<SP(vector_plus(qs,vector_plus(ks,ls)),vector_plus(qs,vector_plus(ks,ls)))*SP(vector_plus(qs,ls),vector_plus(qs,ls))<<endl;
                    cout<<"final result: "<<jac_flow*flow_h*obs_res.jac*integrand_val<<endl;
                };

                obs_res.jac=obs_res.jac*jac_flow*flow_h;

                dreal colour_average=group_average(obs_res.spin1,obs_res.spin2);

                obs_res.eval=colour_average*integrand_val;

                return obs_res;

            };

            observable obs_res_f;
            obs_res_f.eval=0;
            obs_res_f.jac=0;
            obs_res_f.j1={0,0,0,0};
            obs_res_f.j2={0,0,0,0};
            obs_res_f.spin1=0;
            obs_res_f.spin2=0;
            obs_res_f.pg={0,0,0,0};


            return obs_res_f;
            
        };

        observable NLO_t_channel_qq_DrellYan_integrand(vector<dreal> k, vector<dreal> l, vector<dreal> q, dreal a, int comp){
            
            if(comp==1){
                
                vector<vector<dreal>> ni={vector_plus(k,q),vector_swap(vector_plus(k,l))};
                vector<vector<dreal>> nf={vector_swap(l)};

                if(my_flow.t_val(ni,nf,q)==999){
                    observable res;
                    vector<dreal> j1={0,0,0,0};
                    vector<dreal> j2={0,0,0,0};
                    vector<dreal> pg={0,0,0,0};
                    res.eval=0;
                    res.jac=0;
                    res.j1=j1;
                    res.j2=j2;
                    res.spin1=0;
                    res.spin2=0;
                    res.pg=pg;
                    return res;
                };

                vector<dreal> ks=my_flow.perform_flow(ni,nf,q,k);
                vector<dreal> ls=my_flow.perform_flow(ni,nf,q,l);
                vector<dreal> qs=my_flow.perform_flow(ni,nf,q,q);

                dreal jac_flow=my_flow.jacques(ni,nf,q,3);
                dreal flow_h=my_flow.h(my_flow.t_val(ni,nf,q));

                if(DEBUG==1){
                    cout<<"k:  "<<k[0]<<" "<<k[1]<<" "<<k[2]<<endl;
                    cout<<"l:  "<<l[0]<<" "<<l[1]<<" "<<l[2]<<endl;
                    cout<<"q:  "<<q[0]<<" "<<q[1]<<" "<<q[2]<<endl;
                    cout<<"scaled_momenta"<<endl;
                    cout<<"k:  "<<ks[0]<<" "<<ks[1]<<" "<<ks[2]<<endl;
                    cout<<"l:  "<<ls[0]<<" "<<ls[1]<<" "<<ls[2]<<endl;
                    cout<<"q:  "<<qs[0]<<" "<<qs[1]<<" "<<qs[2]<<endl;
                    cout<<"t_value"<<endl;
                    cout<<my_flow.t_val(ni,nf,q)<<endl;
                    cout<<"scaling jacobian"<<endl;
                    cout<<my_flow.jacques(ni,nf,q,3)<<endl;
                };

                
                dreal eq=sqrt(pow(norm(qs),2.)+pow(m,2));
                dreal el=norm(ls);  
                dreal ekl=norm(vector_plus(ks,ls));

                ks.insert(ks.begin(), -ekl+el);
                ls.insert(ls.begin(), -el);
                qs.insert(qs.begin(), eq);

                vector<vector<dreal>> constituents={vector_plus(ks,qs),vector_swap(vector_plus(ks,ls))};

                int spins[2]={-1,1};

                //vector<dreal> obs_res=my_obs.b2b_sampling(constituents, 2, spins, a);
                observable obs_res=my_obs.b2b_sampling_final(constituents, 2, vector_swap(ls), spins, a);

                if(DEBUG==1){
                    vector<dreal> ksp={0.3, 0.1, 0.5};
                    vector<dreal> lsp={0.2, 0.7, 0.1};
                    vector<dreal> qsp={0.1, 0.4, 0.5};

                    dreal ekp=norm(ksp);
                    dreal eqp=sqrt(pow(norm(qsp),2.)+pow(1,2));
                    dreal elp=norm(lsp);  
                    dreal eklp=norm(vector_plus(ksp,lsp));
                    dreal elqp=norm(vector_plus(lsp,qsp));  
                    dreal eklqp=norm(vector_plus(lsp,vector_plus(ksp,qsp)));

                    ksp.insert(ksp.begin(), eklp-elp);
                    lsp.insert(lsp.begin(), elp);
                    qsp.insert(qsp.begin(), eqp);

                    cout<<"m: "<<m<<endl;
                    cout<<"ksp: "<<ksp[0]<<" "<<ksp[1]<<" "<<ksp[2]<<" "<<ksp[3]<<endl;
                    cout<<"lsp: "<<lsp[0]<<" "<<lsp[1]<<" "<<lsp[2]<<" "<<lsp[3]<<endl;
                    cout<<"qsp: "<<qsp[0]<<" "<<qsp[1]<<" "<<qsp[2]<<" "<<qsp[3]<<endl;
                    cout<<"num: "<<NLO_SE_DrellYan_numerator(ksp,vector_swap(vector_plus(lsp,ksp)),vector_swap(qsp))<<endl;
                    cout<<"den: "<<(16*abs(elp*eklp*(-elp+eklp+eqp)*eqp)*SP(ksp,ksp)*SP(ksp,ksp))<<endl;

                    cout<<NLO_SE_DrellYan_numerator(ksp,vector_swap(vector_plus(lsp,ksp)),vector_swap(qsp))/
                    (16*abs(elp*eklp*(-elp+eklp+eqp)*eqp)*SP(ksp,ksp)*SP(ksp,ksp))<<endl;
                };


                dreal integrand_val=constant*NLO_SE_DrellYan_numerator(ks,vector_swap(vector_plus(ls,ks)),vector_swap(qs))/
                    (16*abs(el*ekl*(-el+ekl-eq)*eq)*SP(ks,ks)*SP(ks,ks));

                if(DEBUG==1){
                    cout<<"CHOSEN VALS"<<endl;
                    vector<dreal> qt={0.50415762567199994E+03, 0.11000192113852899E+03, 0.44113190966797782E+03, -0.19789359716207196E+03};
                    vector<dreal> lt={-0.49584237432800018E+03, 0.11000192113852900E+03, 0.44113190966797788E+03, -0.19789359716207190E+03};
                    vector<dreal> kqt={0.50000000000000000E+03, 0.00000000000000000E+00, 0.00000000000000000E+00, 0.50000000000000000E+03};
                    vector<dreal> kt=vector_minus(kqt,qt);

                    vector<dreal> qt1={0.50415762567199999E+03, -0.27957130752081480E+03, -0.38754144149090394E+03, 0.13232979754190944E+03};
                    vector<dreal> lt1={-0.49584237432800018E+03, -0.27957130752081486E+03, -0.38754144149090416E+03, 0.13232979754190944E+03};
                    vector<dreal> kqt1={0.50000000000000000E+03, 0.00000000000000000E+00, 0.00000000000000000E+00, 0.50000000000000000E+03};
                    vector<dreal> kt1=vector_minus(kqt1,qt1);

                    dreal gV=0.1303037631031056;
                    dreal gs=sqrt(0.118*4*M_PI);
                    dreal NC=3;
                    dreal CF=(NC*NC-1)/(2*NC);
                    dreal spin_norm=2;

                    dreal consty=pow(gV,2)*pow(gs,2)*CF*NC/(pow(spin_norm,2.)*pow(NC,2.));

                    dreal val_comp_1=consty*NLO_SE_DrellYan_numerator(kt,vector_swap(vector_plus(lt,kt)),vector_swap(qt))/
                    (SP(kt,kt)*SP(kt,kt));

                    dreal val_comp_2=consty*NLO_SE_DrellYan_numerator(kt1,vector_swap(vector_plus(lt1,kt1)),vector_swap(qt1))/
                    (SP(kt1,kt1)*SP(kt1,kt1));

                    cout<<val_comp_1<<endl;
                    cout<<val_comp_2<<endl;

                };


                if(DEBUG==1){
                    cout<<"############### cut 1 info #############"<<endl;
                    cout<<"original kinematics"<<endl;
                    cout<<"k:  "<<k[0]<<" "<<k[1]<<" "<<k[2]<<endl;
                    cout<<"l:  "<<l[0]<<" "<<l[1]<<" "<<l[2]<<endl;
                    cout<<"q:  "<<q[0]<<" "<<q[1]<<" "<<q[2]<<endl;
                    cout<<"kinematics"<<endl;
                    cout<<"k:  "<<ks[0]<<" "<<ks[1]<<" "<<ks[2]<<" "<<ks[3]<<endl;
                    cout<<"l:  "<<ls[0]<<" "<<ls[1]<<" "<<ls[2]<<" "<<ls[3]<<endl;
                    cout<<"q:  "<<qs[0]<<" "<<qs[1]<<" "<<qs[2]<<" "<<qs[3]<<endl;
                    cout<<"important momenta"<<endl;
                    cout<<"k:  "<<ks[0]+qs[0]<<" "<<ks[1]+qs[1]<<" "<<ks[2]+qs[2]<<" "<<ks[3]+qs[3]<<endl;
                    cout<<"l:  "<<ls[0]<<" "<<ls[1]<<" "<<ls[2]<<" "<<ls[3]<<endl;
                    cout<<"t value: "<<my_flow.t_val(ni,nf,q)<<endl;
                    cout<<"flow jacobian: "<<jac_flow<<endl;
                    cout<<"h function: "<<flow_h<<endl;
                    cout<<"back 2 back jacobian: "<<obs_res.jac<<endl;
                    cout<<"integrand value: "<<integrand_val<<endl;
                    cout<<"product of energies: "<<abs(el*ekl*(-el+ekl+eq)*eq)<<endl;
                    cout<<"product of propagators: "<<SP(ks,ks)*SP(ks,ks)<<endl;
                    cout<<"final result: "<<jac_flow*flow_h*obs_res.jac*integrand_val<<endl;
                };

                obs_res.jac=obs_res.jac*jac_flow*flow_h;

                dreal colour_average=group_average(obs_res.spin1,obs_res.spin2);

                obs_res.eval=colour_average*integrand_val;

                //cout<<"res:    "<<integrand_val*obs_res[1]<<endl;

                return obs_res;

            };

            if(comp==2){
                
                vector<vector<dreal>> ni={vector_plus(k,q),vector_swap(k)};
                vector<vector<dreal>> nf={};

                if(my_flow.t_val(ni,nf,q)==999){
                    observable res;
                    vector<dreal> j1={0,0,0,0};
                    vector<dreal> j2={0,0,0,0};
                    vector<dreal> pg={0,0,0,0};
                    res.eval=0;
                    res.jac=0;
                    res.j1=j1;
                    res.j2=j2;
                    res.spin1=0;
                    res.spin2=0;
                    res.pg=pg;
                    return res;
                };

                vector<dreal> ks=my_flow.perform_flow(ni,nf,q,k);
                vector<dreal> ls=my_flow.perform_flow(ni,nf,q,l);
                vector<dreal> qs=my_flow.perform_flow(ni,nf,q,q);

                dreal jac_flow=my_flow.jacques(ni,nf,q,3);
                dreal flow_h=my_flow.h(my_flow.t_val(ni,nf,q));

                if(DEBUG==1){
                    cout<<"k:  "<<k[0]<<" "<<k[1]<<" "<<k[2]<<endl;
                    cout<<"l:  "<<l[0]<<" "<<l[1]<<" "<<l[2]<<endl;
                    cout<<"q:  "<<q[0]<<" "<<q[1]<<" "<<q[2]<<endl;
                    cout<<"scaled_momenta"<<endl;
                    cout<<"k:  "<<ks[0]<<" "<<ks[1]<<" "<<ks[2]<<endl;
                    cout<<"l:  "<<ls[0]<<" "<<ls[1]<<" "<<ls[2]<<endl;
                    cout<<"q:  "<<qs[0]<<" "<<qs[1]<<" "<<qs[2]<<endl;
                    cout<<"t_value"<<endl;
                    cout<<my_flow.t_val(ni,nf,q)<<endl;
                    cout<<"scaling jacobian"<<endl;
                    cout<<my_flow.jacques(ni,nf,q,3)<<endl;
                };

                dreal ek=norm(ks);
                dreal eq=sqrt(pow(norm(qs),2.)+pow(m,2));
                dreal el=norm(ls);  

                ks.insert(ks.begin(), -ek);
                ls.insert(ls.begin(), -el);
                qs.insert(qs.begin(), eq);

                vector<vector<dreal>> constituents={vector_plus(ks,qs),vector_swap(ks),vector_swap(ls)};
                //vector<vector<dreal>> constituents={vector_plus(ks,qs),vector_swap(ks)};

                int spins[3]={-1,1,0};
                //int spins[3]={-1,1};

                //vector<dreal> obs_res=my_obs.b2b_sampling(constituents, 3, spins, a);
                observable obs_res=my_obs.b2b_sampling_final(constituents, 3, vector_swap(ls), spins, a);
                //vector<dreal> obs_res=my_obs.b2b_sampling_final(constituents, 2, vector_swap(ls), spins, a);

                if(DEBUG==1){

                    vector<dreal> ksp={0.3, 0.1, 0.5};
                    vector<dreal> lsp={0.2, 0.7, 0.1};
                    vector<dreal> qsp={0.1, 0.4, 0.5};

                    dreal ekp=norm(ksp);
                    dreal eqp=sqrt(pow(norm(qsp),2.)+pow(1,2));
                    dreal elp=norm(lsp);  
                    dreal eklp=norm(vector_plus(ksp,lsp));
                    dreal elqp=norm(vector_plus(lsp,qsp));  
                    dreal eklqp=norm(vector_plus(lsp,vector_plus(ksp,qsp)));

                    ksp.insert(ksp.begin(), ekp);
                    lsp.insert(lsp.begin(), elp);
                    qsp.insert(qsp.begin(), eqp);


                    cout<<"m: "<<m<<endl;
                    cout<<"ksp: "<<ksp[0]<<" "<<ksp[1]<<" "<<ksp[2]<<" "<<ksp[3]<<endl;
                    cout<<"lsp: "<<lsp[0]<<" "<<lsp[1]<<" "<<lsp[2]<<" "<<lsp[3]<<endl;
                    cout<<"qsp: "<<qsp[0]<<" "<<qsp[1]<<" "<<qsp[2]<<" "<<qsp[3]<<endl;
                    cout<<"num: "<<NLO_SE_DrellYan_numerator(ksp,vector_swap(vector_plus(lsp,ksp)),vector_swap(qsp))<<endl;
                    cout<<"num_d: "<<NLO_SE_DrellYan_numerator_d(ksp,vector_swap(vector_plus(lsp,ksp)),vector_swap(qsp))<<endl;
                    cout<<"energies :"<<32*abs(ekp*ekp*elp*(ekp+eqp)*eqp)<<endl;

                    cout<<NLO_SE_DrellYan_numerator_d(ksp,vector_swap(vector_plus(lsp,ksp)),vector_swap(qsp))/
                    (32*abs(ekp*ekp*elp*(ekp+eqp)*eqp)*SP(vector_plus(ksp,lsp),vector_plus(ksp,lsp)))+
                    NLO_SE_DrellYan_numerator(ksp,vector_swap(vector_plus(lsp,ksp)),vector_swap(qsp))*2*(elp+ekp)/
                    (32*abs(ekp*ekp*elp*(ekp+eqp)*eqp)*SP(vector_plus(ksp,lsp),vector_plus(ksp,lsp))*SP(vector_plus(ksp,lsp),vector_plus(ksp,lsp)))<<endl;
                };

                //SE DERIVATIVE ROUTING

                /*dreal integrand_val=-constant*(NLO_SE_DrellYan_numerator_d(ks,vector_swap(vector_plus(ls,ks)),vector_swap(qs))/
                    (32*abs(ek*ek*el*(-ek+eq)*eq)*SP(vector_plus(ks,ls),vector_plus(ks,ls)))+
                    NLO_SE_DrellYan_numerator(ks,vector_swap(vector_plus(ls,ks)),vector_swap(qs))*2*(el+ek)/
                    (32*abs(ek*ek*el*(-ek+eq)*eq)*SP(vector_plus(ks,ls),vector_plus(ks,ls))*SP(vector_plus(ks,ls),vector_plus(ks,ls))));*/

                dreal integrand_val=-constant*(-NLO_SE_DrellYan_numerator_d(ks,vector_swap(vector_plus(ls,ks)),vector_swap(qs))/
                    (32*abs(ek*ek*el*(-ek+eq)*eq)*SP(vector_plus(ks,ls),vector_plus(ks,ls)))+
                    NLO_SE_DrellYan_numerator(ks,vector_swap(vector_plus(ls,ks)),vector_swap(qs))*2*(el+ek)/
                    (32*abs(ek*ek*el*(-ek+eq)*eq)*SP(vector_plus(ks,ls),vector_plus(ks,ls))*SP(vector_plus(ks,ls),vector_plus(ks,ls))));


                /*cout<<NLO_SE_DrellYan_numerator_d(ks,vector_swap(vector_plus(ls,ks)),vector_swap(qs))/
                    (32*abs(ek*ek*el*(-ek+eq)*eq)*SP(vector_plus(ks,ls),vector_plus(ks,ls)))<<endl;
                cout<<NLO_SE_DrellYan_numerator(ks,vector_swap(vector_plus(ls,ks)),vector_swap(qs))*2*(el+ek)/
                    (32*abs(ek*ek*el*(-ek+eq)*eq)*SP(vector_plus(ks,ls),vector_plus(ks,ls))*SP(vector_plus(ks,ls),vector_plus(ks,ls)))<<endl;
                cout<<(NLO_SE_DrellYan_numerator_d(ks,vector_swap(vector_plus(ls,ks)),vector_swap(qs))/
                    (32*abs(ek*ek*el*(-ek+eq)*eq)*SP(vector_plus(ks,ls),vector_plus(ks,ls)))-
                    NLO_SE_DrellYan_numerator(ks,vector_swap(vector_plus(ls,ks)),vector_swap(qs))*2*(el+ek)/
                    (32*abs(ek*ek*el*(-ek+eq)*eq)*SP(vector_plus(ks,ls),vector_plus(ks,ls))*SP(vector_plus(ks,ls),vector_plus(ks,ls))))<<endl;
                cout<<"k "<<k[0]<<" "<<k[1]<<" "<<k[2]<<endl;
                cout<<"l "<<l[0]<<" "<<l[1]<<" "<<l[2]<<endl;
                cout<<"q "<<q[0]<<" "<<q[1]<<" "<<q[2]<<endl;
                cout<<"t "<<my_flow.t_val(ni,nf,q)<<endl;
                cout<<"Es "<<1/abs(32*ek*ek*el*(-ek+eq)*eq)<<endl;
                cout<<"ek "<<ek<<endl;
                cout<<"el "<<el<<endl;
                cout<<"eq "<<eq<<endl;*/
                
                
                if(DEBUG==1){
                    cout<<"############### cut 2 info #############"<<endl;
                    cout<<"original kinematics"<<endl;
                    cout<<"k:  "<<k[0]<<" "<<k[1]<<" "<<k[2]<<endl;
                    cout<<"l:  "<<l[0]<<" "<<l[1]<<" "<<l[2]<<endl;
                    cout<<"q:  "<<q[0]<<" "<<q[1]<<" "<<q[2]<<endl;
                    cout<<"kinematics"<<endl;
                    cout<<"k:  "<<ks[0]<<" "<<ks[1]<<" "<<ks[2]<<" "<<ks[3]<<endl;
                    cout<<"l:  "<<ls[0]<<" "<<ls[1]<<" "<<ls[2]<<" "<<ls[3]<<endl;
                    cout<<"q:  "<<qs[0]<<" "<<qs[1]<<" "<<qs[2]<<" "<<qs[3]<<endl;
                    cout<<"important momenta"<<endl;
                    cout<<"k:  "<<ks[0]+qs[0]<<" "<<ks[1]+qs[1]<<" "<<ks[2]+qs[2]<<" "<<ks[3]+qs[3]<<endl;
                    cout<<"l:  "<<ls[0]<<" "<<ls[1]<<" "<<ls[2]<<" "<<ls[3]<<endl;
                    cout<<"t value: "<<my_flow.t_val(ni,nf,q)<<endl;
                    cout<<"flow jacobian: "<<jac_flow<<endl;
                    cout<<"h function: "<<flow_h<<endl;
                    cout<<"back 2 back jacobian: "<<obs_res.jac<<endl;
                    cout<<"integrand value: "<<integrand_val<<endl;
                    cout<<"product of energies: "<<abs(ek*ek*el*(ek+eq)*eq)<<endl;
                    cout<<"product of propagators: "<<SP(vector_plus(ks,ls),vector_plus(ks,ls))<<endl;
                    cout<<"final result: "<<jac_flow*flow_h*obs_res.jac*integrand_val<<endl;
                };


                obs_res.jac=obs_res.jac*jac_flow*flow_h;

                dreal colour_average=group_average(obs_res.spin1,obs_res.spin2);

                obs_res.eval=colour_average*integrand_val;

                return obs_res;

            };

            observable obs_res_f;
            obs_res_f.eval=0;
            obs_res_f.jac=0;
            obs_res_f.j1={0,0,0,0};
            obs_res_f.j2={0,0,0,0};
            obs_res_f.spin1=0;
            obs_res_f.spin2=0;
            obs_res_f.pg={0,0,0,0};


            return obs_res_f;
            
        };

        observable NLO_t_channel_qg_DrellYan_integrand(vector<dreal> k, vector<dreal> l, vector<dreal> q, dreal a, int comp){
            
            if(comp==1){
                
                vector<vector<dreal>> ni={vector_plus(k,q),vector_swap(vector_plus(k,l))};
                vector<vector<dreal>> nf={vector_swap(l)};

                if(my_flow.t_val(ni,nf,q)==999){
                    observable res;
                    vector<dreal> j1={0,0,0,0};
                    vector<dreal> j2={0,0,0,0};
                    vector<dreal> pg={0,0,0,0};
                    res.eval=0;
                    res.jac=0;
                    res.j1=j1;
                    res.j2=j2;
                    res.spin1=0;
                    res.spin2=0;
                    res.pg=pg;
                    return res;
                };

                vector<dreal> ks=my_flow.perform_flow(ni,nf,q,k);
                vector<dreal> ls=my_flow.perform_flow(ni,nf,q,l);
                vector<dreal> qs=my_flow.perform_flow(ni,nf,q,q);

                dreal jac_flow=my_flow.jacques(ni,nf,q,3);
                dreal flow_h=my_flow.h(my_flow.t_val(ni,nf,q));

                if(DEBUG==1){
                    cout<<"k:  "<<k[0]<<" "<<k[1]<<" "<<k[2]<<endl;
                    cout<<"l:  "<<l[0]<<" "<<l[1]<<" "<<l[2]<<endl;
                    cout<<"q:  "<<q[0]<<" "<<q[1]<<" "<<q[2]<<endl;
                    cout<<"scaled_momenta"<<endl;
                    cout<<"k:  "<<ks[0]<<" "<<ks[1]<<" "<<ks[2]<<endl;
                    cout<<"l:  "<<ls[0]<<" "<<ls[1]<<" "<<ls[2]<<endl;
                    cout<<"q:  "<<qs[0]<<" "<<qs[1]<<" "<<qs[2]<<endl;
                    cout<<"t_value"<<endl;
                    cout<<my_flow.t_val(ni,nf,q)<<endl;
                    cout<<"scaling jacobian"<<endl;
                    cout<<my_flow.jacques(ni,nf,q,3)<<endl;
                };

                
                dreal eq=sqrt(pow(norm(qs),2.)+pow(m,2));
                dreal el=norm(ls);  
                dreal ekl=norm(vector_plus(ks,ls));

                ks.insert(ks.begin(), -ekl+el);
                ls.insert(ls.begin(), -el);
                qs.insert(qs.begin(), eq);

                vector<vector<dreal>> constituents={vector_plus(ks,qs),vector_swap(vector_plus(ks,ls))};

                int spins[2]={-1,0};

                observable obs_res=my_obs.b2b_sampling_final(constituents, 2, vector_swap(ls), spins, a);

                if(DEBUG==1){
                    vector<dreal> ksp={0.3, 0.1, 0.5};
                    vector<dreal> lsp={0.2, 0.7, 0.1};
                    vector<dreal> qsp={0.1, 0.4, 0.5};

                    dreal ekp=norm(ksp);
                    dreal eqp=sqrt(pow(norm(qsp),2.)+pow(m,2));
                    dreal elp=norm(lsp);  
                    dreal eklp=norm(vector_plus(ksp,lsp));
                    dreal elqp=norm(vector_plus(lsp,qsp));  
                    dreal eklqp=norm(vector_plus(lsp,vector_plus(ksp,qsp)));

                    ksp.insert(ksp.begin(), -eklp+elp);
                    lsp.insert(lsp.begin(), -elp);
                    qsp.insert(qsp.begin(), eqp);


                    cout<<"HEYYYY"<<endl;
                    cout<<"m: "<<m<<endl;
                    cout<<"ksp: "<<ksp[0]<<" "<<ksp[1]<<" "<<ksp[2]<<" "<<ksp[3]<<endl;
                    cout<<"lsp: "<<lsp[0]<<" "<<lsp[1]<<" "<<lsp[2]<<" "<<lsp[3]<<endl;
                    cout<<"qsp: "<<qsp[0]<<" "<<qsp[1]<<" "<<qsp[2]<<" "<<qsp[3]<<endl;
                    cout<<"num: "<<NLO_SE_DrellYan_numerator(vector_swap(ksp),lsp,qsp)<<endl;
                    cout<<"den: "<<(16*abs(elp*eklp*(-eklp+elp+eqp)*eqp)*SP(ksp,ksp)*SP(ksp,ksp))<<endl;

                    cout<<NLO_SE_DrellYan_numerator(vector_swap(ksp),lsp,qsp)/
                    (16*abs(elp*eklp*(-eklp+elp+eqp)*eqp)*SP(ksp,ksp)*SP(ksp,ksp))<<endl;
                };

                dreal integrand_val=constant*NLO_SE_DrellYan_numerator(vector_swap(ks),ls,qs)/
                    (16*abs(el*ekl*(-ekl+el+eq)*eq)*SP(ks,ks)*SP(ks,ks));

                //cout<<integrand_val<<endl;
                
                if(DEBUG==1){  
                    vector<dreal> lsp={-104.621, 22.4313, 15.3108, 101.035};
                    vector<dreal> qsp={138.784, 22.4313, 15.3108, 101.035};
                    vector<dreal> kqsp={-121.702, 0., 0., -121.702};
                    vector<dreal> ksp=vector_minus(kqsp,qsp);

                    cout<<"comparison_MG: "<<constant*NLO_SE_DrellYan_numerator(vector_swap(ksp),lsp,qsp)/(SP(ksp,ksp)*SP(ksp,ksp))<<endl;
                };

                if(DEBUG==1){
                    vector<dreal> lsp={-41.6672, -37.4207, 18.0886, 2.94203};
                    vector<dreal> qsp={100.25671940723308, -37.4207, 18.0886, 2.94203};
                    vector<dreal> kqsp={-70.962, 0., 0., -70.962};
                    vector<dreal> ksp=vector_minus(kqsp,qsp);

                    cout<<"comparison_MG_2: "<<constant*NLO_SE_DrellYan_numerator(vector_swap(ksp),lsp,qsp)/(SP(ksp,ksp)*SP(ksp,ksp))<<endl;
                };

                if(DEBUG==1){
                    cout<<"############### cut 1 info #############"<<endl;
                    cout<<"original kinematics"<<endl;
                    cout<<"k:  "<<k[0]<<" "<<k[1]<<" "<<k[2]<<endl;
                    cout<<"l:  "<<l[0]<<" "<<l[1]<<" "<<l[2]<<endl;
                    cout<<"q:  "<<q[0]<<" "<<q[1]<<" "<<q[2]<<endl;
                    cout<<"ekl: "<<ekl<<endl;
                    cout<<"kinematics"<<endl;
                    cout<<"k:  "<<ks[0]<<" "<<ks[1]<<" "<<ks[2]<<" "<<ks[3]<<endl;
                    cout<<"l:  "<<ls[0]<<" "<<ls[1]<<" "<<ls[2]<<" "<<ls[3]<<endl;
                    cout<<"q:  "<<qs[0]<<" "<<qs[1]<<" "<<qs[2]<<" "<<qs[3]<<endl;
                    cout<<"important momenta"<<endl;
                    cout<<"k:  "<<ks[0]+qs[0]<<" "<<ks[1]+qs[1]<<" "<<ks[2]+qs[2]<<" "<<ks[3]+qs[3]<<endl;
                    cout<<"l:  "<<ls[0]<<" "<<ls[1]<<" "<<ls[2]<<" "<<ls[3]<<endl;
                    cout<<"t value: "<<my_flow.t_val(ni,nf,q)<<endl;
                    cout<<"flow jacobian: "<<jac_flow<<endl;
                    cout<<"h function: "<<flow_h<<endl;
                    cout<<"back 2 back jacobian: "<<obs_res.jac<<endl;
                    cout<<"integrand value: "<<integrand_val<<endl;
                    cout<<"product of energies: "<<abs(el*ekl*(-ekl+el+eq)*eq)<<endl;
                    cout<<"product of propagators: "<<SP(ks,ks)*SP(ks,ks)<<endl;
                    cout<<"final result: "<<jac_flow*flow_h*obs_res.jac*integrand_val<<endl;
                };

                obs_res.jac=obs_res.jac*jac_flow*flow_h;

                dreal colour_average=group_average(obs_res.spin1,obs_res.spin2);


                obs_res.eval=colour_average*integrand_val;

                return obs_res;

            };

            if(comp==2){
                
                vector<vector<dreal>> ni={vector_plus(k,q),vector_swap(k)};
                vector<vector<dreal>> nf={};

                if(my_flow.t_val(ni,nf,q)==999){
                    observable res;
                    vector<dreal> j1={0,0,0,0};
                    vector<dreal> j2={0,0,0,0};
                    vector<dreal> pg={0,0,0,0};
                    res.eval=0;
                    res.jac=0;
                    res.j1=j1;
                    res.j2=j2;
                    res.spin1=0;
                    res.spin2=0;
                    res.pg=pg;
                    return res;
                };

                vector<dreal> ks=my_flow.perform_flow(ni,nf,q,k);
                vector<dreal> ls=my_flow.perform_flow(ni,nf,q,l);
                vector<dreal> qs=my_flow.perform_flow(ni,nf,q,q);

                dreal jac_flow=my_flow.jacques(ni,nf,q,3);
                dreal flow_h=my_flow.h(my_flow.t_val(ni,nf,q));

                if(DEBUG==1){
                    cout<<"k:  "<<k[0]<<" "<<k[1]<<" "<<k[2]<<endl;
                    cout<<"l:  "<<l[0]<<" "<<l[1]<<" "<<l[2]<<endl;
                    cout<<"q:  "<<q[0]<<" "<<q[1]<<" "<<q[2]<<endl;
                    cout<<"scaled_momenta"<<endl;
                    cout<<"k:  "<<ks[0]<<" "<<ks[1]<<" "<<ks[2]<<endl;
                    cout<<"l:  "<<ls[0]<<" "<<ls[1]<<" "<<ls[2]<<endl;
                    cout<<"q:  "<<qs[0]<<" "<<qs[1]<<" "<<qs[2]<<endl;
                    cout<<"t_value"<<endl;
                    cout<<my_flow.t_val(ni,nf,q)<<endl;
                    cout<<"scaling jacobian"<<endl;
                    cout<<my_flow.jacques(ni,nf,q,3)<<endl;
                };

                dreal ek=norm(ks);
                dreal eq=sqrt(pow(norm(qs),2.)+pow(m,2));
                dreal el=norm(ls);  

                ks.insert(ks.begin(), -ek);
                ls.insert(ls.begin(), -el);
                qs.insert(qs.begin(), eq);

                vector<vector<dreal>> constituents={vector_plus(ks,qs),vector_swap(ks),vector_swap(ls)};

                int spins[3]={-1,1,-1};

                observable obs_res=my_obs.b2b_sampling_final(constituents, 3, vector_swap(ls), spins, a);

                if(DEBUG==1){
                    vector<dreal> ksp={0.3, 0.1, 0.5};
                    vector<dreal> lsp={0.2, 0.7, 0.1};
                    vector<dreal> qsp={0.1, 0.4, 0.5};

                    dreal ekp=norm(ksp);
                    dreal eqp=sqrt(pow(norm(qsp),2.)+pow(m,2));
                    dreal elp=norm(lsp);  
                    dreal eklp=norm(vector_plus(ksp,lsp));
                    dreal elqp=norm(vector_plus(lsp,qsp));  
                    dreal eklqp=norm(vector_plus(lsp,vector_plus(ksp,qsp)));

                    ksp.insert(ksp.begin(), -ekp);
                    lsp.insert(lsp.begin(), -elp);
                    qsp.insert(qsp.begin(), eqp);


                    cout<<"m: "<<m<<endl;
                    cout<<"ksp: "<<ksp[0]<<" "<<ksp[1]<<" "<<ksp[2]<<" "<<ksp[3]<<endl;
                    cout<<"lsp: "<<lsp[0]<<" "<<lsp[1]<<" "<<lsp[2]<<" "<<lsp[3]<<endl;
                    cout<<"qsp: "<<qsp[0]<<" "<<qsp[1]<<" "<<qsp[2]<<" "<<qsp[3]<<endl;
                    cout<<"num: "<<NLO_SE_DrellYan_numerator_d(vector_swap(ksp),lsp,qsp)<<endl;
                    cout<<"num_d: "<<NLO_SE_DrellYan_numerator(vector_swap(ksp),lsp,qsp)*2*(ekp+elp)<<endl;
                    cout<<"energies :"<<32*abs(ekp*ekp*elp*(-ekp+eqp)*eqp)<<endl;

                    cout<</*NLO_SE_DrellYan_numerator_d(vector_swap(ksp),lsp,qsp)/
                    (32*abs(ekp*ekp*elp*(-ekp+eqp)*eqp)*SP(vector_plus(ksp,lsp),vector_plus(ksp,lsp)))-*/
                    NLO_SE_DrellYan_numerator(vector_swap(ksp),lsp,qsp)*2*(elp+ekp)/
                    (32*abs(ekp*ekp*elp*(-ekp+eqp)*eqp)*SP(vector_plus(ksp,lsp),vector_plus(ksp,lsp))*SP(vector_plus(ksp,lsp),vector_plus(ksp,lsp)))<<endl;
                };


                dreal integrand_val=/*NLO_SE_DrellYan_numerator_d(vector_swap(ks),ls,qs)/
                    (32*abs(ek*ek*el*(-ek+eq)*eq)*SP(vector_plus(ks,ls),vector_plus(ks,ls)))+*/
                    -constant*NLO_SE_DrellYan_numerator(vector_swap(ks),ls,qs)*2*(el+ek)/
                    (32*abs(ek*ek*el*(-ek+eq)*eq)*SP(vector_plus(ks,ls),vector_plus(ks,ls))*SP(vector_plus(ks,ls),vector_plus(ks,ls)));

                //cout<<integrand_val<<endl;

                
                if(DEBUG==1){
                    cout<<"############### cut 2 info #############"<<endl;
                    cout<<"original kinematics"<<endl;
                    cout<<"k:  "<<k[0]<<" "<<k[1]<<" "<<k[2]<<endl;
                    cout<<"l:  "<<l[0]<<" "<<l[1]<<" "<<l[2]<<endl;
                    cout<<"q:  "<<q[0]<<" "<<q[1]<<" "<<q[2]<<endl;
                    cout<<"kinematics"<<endl;
                    cout<<"k:  "<<ks[0]<<" "<<ks[1]<<" "<<ks[2]<<" "<<ks[3]<<endl;
                    cout<<"l:  "<<ls[0]<<" "<<ls[1]<<" "<<ls[2]<<" "<<ls[3]<<endl;
                    cout<<"q:  "<<qs[0]<<" "<<qs[1]<<" "<<qs[2]<<" "<<qs[3]<<endl;
                    cout<<"important momenta"<<endl;
                    cout<<"k+l:  "<<ks[0]+ls[0]<<" "<<ks[1]+ls[1]<<" "<<ks[2]+ls[2]<<" "<<ks[3]+ls[3]<<endl;
                    cout<<"l:  "<<ls[0]<<" "<<ls[1]<<" "<<ls[2]<<" "<<ls[3]<<endl;
                    cout<<"t value: "<<my_flow.t_val(ni,nf,q)<<endl;
                    cout<<"flow jacobian: "<<jac_flow<<endl;
                    cout<<"h function: "<<flow_h<<endl;
                    cout<<"back 2 back jacobian: "<<obs_res.jac<<endl;
                    cout<<"integrand value: "<<integrand_val<<endl;
                    cout<<"product of energies: "<<abs(ek*ek*el*(ek+eq)*eq)<<endl;
                    cout<<"ek,el,ekl: "<<ek<<" "<<el<<" "<<norm(vector_plus(ks,ls))<<"    surf: "<<ek+el-norm(vector_plus(ks,ls))<<endl;
                    cout<<"product of propagators: "<<SP(vector_plus(ks,ls),vector_plus(ks,ls))*SP(vector_plus(ks,ls),vector_plus(ks,ls))<<endl;
                    cout<<"final result: "<<jac_flow*flow_h*obs_res.jac*integrand_val<<endl;
                };

                obs_res.jac=obs_res.jac*jac_flow*flow_h;

                dreal colour_average=group_average(obs_res.spin1,obs_res.spin2);


                obs_res.eval=colour_average*integrand_val;

                return obs_res;

            };

            observable obs_res_f;
            obs_res_f.eval=0;
            obs_res_f.jac=0;
            obs_res_f.j1={0,0,0,0};
            obs_res_f.j2={0,0,0,0};
            obs_res_f.spin1=0;
            obs_res_f.spin2=0;
            obs_res_f.pg={0,0,0,0};


            return obs_res_f;
            
        };

        observable NLO_u_channel_DrellYan_integrand(vector<dreal> k, vector<dreal> l, vector<dreal> q, dreal a, int comp){
            

            if(comp==1){

                if(DEBUG==1){cout<<"-------------------COMPONENT: "<<comp<<" -------------------------"<<endl;};
                
                vector<vector<dreal>> ni={vector_plus(k,q),vector_swap(vector_plus(k,l))};
                vector<vector<dreal>> nf={vector_swap(l)};

                

                if(my_flow.t_val(ni,nf,q)==999){
                    observable res;
                    vector<dreal> j1={0,0,0,0};
                    vector<dreal> j2={0,0,0,0};
                    vector<dreal> pg={0,0,0,0};
                    res.eval=0;
                    res.jac=0;
                    res.j1=j1;
                    res.j2=j2;
                    res.spin1=0;
                    res.spin2=0;
                    res.pg=pg;
                    return res;
                };

                vector<dreal> ks=my_flow.perform_flow(ni,nf,q,k);
                vector<dreal> ls=my_flow.perform_flow(ni,nf,q,l);
                vector<dreal> qs=my_flow.perform_flow(ni,nf,q,q);

                dreal jac_flow=my_flow.jacques(ni,nf,q,3);
                dreal flow_h=my_flow.h(my_flow.t_val(ni,nf,q));

                if(DEBUG==1){
                    cout<<"k:  "<<k[0]<<" "<<k[1]<<" "<<k[2]<<endl;
                    cout<<"l:  "<<l[0]<<" "<<l[1]<<" "<<l[2]<<endl;
                    cout<<"q:  "<<q[0]<<" "<<q[1]<<" "<<q[2]<<endl;
                    cout<<"scaled_momenta"<<endl;
                    cout<<"k:  "<<ks[0]<<" "<<ks[1]<<" "<<ks[2]<<endl;
                    cout<<"l:  "<<ls[0]<<" "<<ls[1]<<" "<<ls[2]<<endl;
                    cout<<"q:  "<<qs[0]<<" "<<qs[1]<<" "<<qs[2]<<endl;
                    cout<<"t_value"<<endl;
                    cout<<my_flow.t_val(ni,nf,q)<<endl;
                    cout<<"scaling jacobian"<<endl;
                    cout<<my_flow.jacques(ni,nf,q,3)<<endl;
                };

                
                dreal eq=sqrt(pow(norm(qs),2.)+pow(m,2));
                dreal el=norm(ls);  
                dreal ekl=norm(vector_plus(ks,ls));

                ks.insert(ks.begin(), -ekl+el);
                ls.insert(ls.begin(), -el);
                qs.insert(qs.begin(), eq);

                vector<vector<dreal>> constituents={vector_plus(ks,qs),vector_swap(vector_plus(ks,ls))};

                int spins[2]={-1,1};

                observable obs_res=my_obs.b2b_sampling_final(constituents, 2, vector_swap(ls), spins, a);

                if(DEBUG==1){
                    vector<dreal> ksp={0.3, 0.1, 0.5};
                    vector<dreal> lsp={0.2, 0.7, 0.1};
                    vector<dreal> qsp={0.1, 0.4, 0.5};

                    dreal ekp=norm(ksp);
                    dreal eqp=sqrt(pow(norm(qsp),2.)+pow(m,2));
                    dreal elp=norm(lsp);  
                    dreal eklp=norm(vector_plus(ksp,lsp));
                    dreal elqp=norm(vector_plus(lsp,qsp));  
                    dreal eklqp=norm(vector_plus(lsp,vector_plus(ksp,qsp)));

                    ksp.insert(ksp.begin(), -eklp+elp);
                    lsp.insert(lsp.begin(), -elp);
                    qsp.insert(qsp.begin(), eqp);


                    //NUMERATOR HERE IS OK

                    cout<<"m: "<<m<<endl;
                    cout<<"ksp: "<<ksp[0]<<" "<<ksp[1]<<" "<<ksp[2]<<" "<<ksp[3]<<endl;
                    cout<<"lsp: "<<lsp[0]<<" "<<lsp[1]<<" "<<lsp[2]<<" "<<lsp[3]<<endl;
                    cout<<"qsp: "<<qsp[0]<<" "<<qsp[1]<<" "<<qsp[2]<<" "<<qsp[3]<<endl;
                    cout<<"num: "<<NLO_DT_DrellYan_numerator(vector_swap(vector_plus(ksp,lsp)),vector_swap(ksp),qsp)<<endl;
                    cout<<"den: "<<(16*abs(elp*eklp*(eklp-elp+eqp)*eqp)*SP(ksp,ksp)*SP(vector_plus(ksp,vector_plus(lsp,qsp)),vector_plus(ksp,vector_plus(lsp,qsp))))<<endl;

                    cout<<NLO_DT_DrellYan_numerator(vector_swap(vector_plus(ksp,lsp)),vector_swap(ksp),qsp)/
                    (16*abs(elp*eklp*(eklp-elp+eqp)*eqp)*SP(ksp,ksp)*SP(vector_plus(ksp,vector_plus(lsp,qsp)),vector_plus(ksp,vector_plus(lsp,qsp))))<<endl;
                };


                dreal integrand_val=2*constant*NLO_DT_DrellYan_numerator(vector_swap(vector_plus(ks,ls)),vector_swap(ks),qs)/
                    (16*abs(el*ekl*(-ekl+el+eq)*eq)*SP(ks,ks)*SP(vector_plus(ks,vector_plus(ls,qs)),vector_plus(ks,vector_plus(ls,qs))));
                
                if(DEBUG==1){
                    cout<<"############### cut 1 info #############"<<endl;
                    cout<<"original kinematics"<<endl;
                    cout<<"k:  "<<k[0]<<" "<<k[1]<<" "<<k[2]<<endl;
                    cout<<"l:  "<<l[0]<<" "<<l[1]<<" "<<l[2]<<endl;
                    cout<<"q:  "<<q[0]<<" "<<q[1]<<" "<<q[2]<<endl;
                    cout<<"kinematics"<<endl;
                    cout<<"k:  "<<ks[0]<<" "<<ks[1]<<" "<<ks[2]<<" "<<ks[3]<<endl;
                    cout<<"l:  "<<ls[0]<<" "<<ls[1]<<" "<<ls[2]<<" "<<ls[3]<<endl;
                    cout<<"q:  "<<qs[0]<<" "<<qs[1]<<" "<<qs[2]<<" "<<qs[3]<<endl;
                    cout<<"important momenta"<<endl;
                    cout<<"k:  "<<ks[0]+qs[0]<<" "<<ks[1]+qs[1]<<" "<<ks[2]+qs[2]<<" "<<ks[3]+qs[3]<<endl;
                    cout<<"l:  "<<ls[0]<<" "<<ls[1]<<" "<<ls[2]<<" "<<ls[3]<<endl;
                    cout<<"flow check: "<<-ekl+el+eq<<" "<<my_flow.t_val(ni,nf,q)*norm(vector_plus(k,q))<<endl;
                    cout<<"t value: "<<my_flow.t_val(ni,nf,q)<<endl;
                    cout<<"flow jacobian: "<<jac_flow<<endl;
                    cout<<"h function: "<<flow_h<<endl;
                    cout<<"back 2 back jacobian: "<<obs_res.jac<<endl;
                    cout<<"integrand value: "<<integrand_val<<endl;
                    cout<<"product of energies: "<<abs(el*ekl*(-ekl+el+eq)*eq)<<endl;
                    cout<<"product of propagators: "<<SP(ks,ks)*SP(vector_plus(ks,vector_plus(ls,qs)),vector_plus(ks,vector_plus(ls,qs)))<<endl;
                    cout<<"final result: "<<jac_flow*flow_h*obs_res.jac*integrand_val<<endl;
                };

                obs_res.jac=obs_res.jac*jac_flow*flow_h;

                dreal colour_average=group_average(obs_res.spin1,obs_res.spin2);

                obs_res.eval=colour_average*integrand_val;

                return obs_res;

            };

            if(comp==2){

                if(DEBUG==1){cout<<"-------------------COMPONENT: "<<comp<<" -------------------------"<<endl;};
                
                vector<vector<dreal>> ni={vector_plus(k,q),vector_swap(k)};
                vector<vector<dreal>> nf={};

                if(my_flow.t_val(ni,nf,q)==999){
                    observable res;
                    vector<dreal> j1={0,0,0,0};
                    vector<dreal> j2={0,0,0,0};
                    vector<dreal> pg={0,0,0,0};
                    res.eval=0;
                    res.jac=0;
                    res.j1=j1;
                    res.j2=j2;
                    res.spin1=0;
                    res.spin2=0;
                    res.pg=pg;
                    return res;
                };

                vector<dreal> ks=my_flow.perform_flow(ni,nf,q,k);
                vector<dreal> ls=my_flow.perform_flow(ni,nf,q,l);
                vector<dreal> qs=my_flow.perform_flow(ni,nf,q,q);

                dreal jac_flow=my_flow.jacques(ni,nf,q,3);
                dreal flow_h=my_flow.h(my_flow.t_val(ni,nf,q));

                if(DEBUG==1){
                    cout<<"k:  "<<k[0]<<" "<<k[1]<<" "<<k[2]<<endl;
                    cout<<"l:  "<<l[0]<<" "<<l[1]<<" "<<l[2]<<endl;
                    cout<<"q:  "<<q[0]<<" "<<q[1]<<" "<<q[2]<<endl;
                    cout<<"scaled_momenta"<<endl;
                    cout<<"k:  "<<ks[0]<<" "<<ks[1]<<" "<<ks[2]<<endl;
                    cout<<"l:  "<<ls[0]<<" "<<ls[1]<<" "<<ls[2]<<endl;
                    cout<<"q:  "<<qs[0]<<" "<<qs[1]<<" "<<qs[2]<<endl;
                    cout<<"t_value"<<endl;
                    cout<<my_flow.t_val(ni,nf,q)<<endl;
                    cout<<"scaling jacobian"<<endl;
                    cout<<my_flow.jacques(ni,nf,q,3)<<endl;
                };

                
                dreal eq=sqrt(pow(norm(qs),2.)+pow(m,2));
                dreal el=norm(ls);  
                dreal ek=norm(ks);

                ks.insert(ks.begin(), -ek);
                ls.insert(ls.begin(), -el);
                qs.insert(qs.begin(), eq);

                vector<vector<dreal>> constituents={vector_plus(ks,qs),vector_swap(ks),vector_swap(ls)};

                int spins[3]={-1,1,0};

                observable obs_res=my_obs.b2b_sampling_final(constituents, 3, vector_swap(ls), spins, a);

                if(DEBUG==1){
                    vector<dreal> ksp={0.3, 0.1, 0.5};
                    vector<dreal> lsp={0.2, 0.7, 0.1};
                    vector<dreal> qsp={0.1, 0.4, 0.5};

                    dreal ekp=norm(ksp);
                    dreal eqp=sqrt(pow(norm(qsp),2.)+pow(m,2));
                    dreal elp=norm(lsp);  
                    dreal eklp=norm(vector_plus(ksp,lsp));
                    dreal elqp=norm(vector_plus(lsp,qsp));  
                    dreal eklqp=norm(vector_plus(lsp,vector_plus(ksp,qsp)));

                    ksp.insert(ksp.begin(), -ekp);
                    lsp.insert(lsp.begin(), -elp);
                    qsp.insert(qsp.begin(), eqp);

                    //NUMERATOR HERE IS OK

                    cout<<"m: "<<m<<endl;
                    cout<<"ksp: "<<ksp[0]<<" "<<ksp[1]<<" "<<ksp[2]<<" "<<ksp[3]<<endl;
                    cout<<"lsp: "<<lsp[0]<<" "<<lsp[1]<<" "<<lsp[2]<<" "<<lsp[3]<<endl;
                    cout<<"qsp: "<<qsp[0]<<" "<<qsp[1]<<" "<<qsp[2]<<" "<<qsp[3]<<endl;
                    cout<<"num: "<<NLO_DT_DrellYan_numerator(vector_swap(vector_plus(ksp,lsp)),vector_swap(ksp),qsp)<<endl;
                    cout<<"den: "<<(16*abs(elp*ekp*(-ekp+eqp)*eqp)*SP(vector_plus(ksp,lsp),vector_plus(ksp,lsp))*SP(vector_plus(ksp,vector_plus(lsp,qsp)),vector_plus(ksp,vector_plus(lsp,qsp))))<<endl;

                    cout<<NLO_DT_DrellYan_numerator(vector_swap(vector_plus(ksp,lsp)),vector_swap(ksp),qsp)/
                    (16*abs(elp*ekp*(-ekp+eqp)*eqp)*SP(vector_plus(ksp,lsp),vector_plus(ksp,lsp))*SP(vector_plus(ksp,vector_plus(lsp,qsp)),vector_plus(ksp,vector_plus(lsp,qsp))))<<endl;
                };


                dreal integrand_val=2*constant*NLO_DT_DrellYan_numerator(vector_swap(vector_plus(ks,ls)),vector_swap(ks),qs)/
                    (16*abs(el*ek*(-ek+eq)*eq)*SP(vector_plus(ks,ls),vector_plus(ks,ls))*SP(vector_plus(ks,vector_plus(ls,qs)),vector_plus(ks,vector_plus(ls,qs))));
                
                if(DEBUG==1){
                    cout<<"############### cut 1 info #############"<<endl;
                    cout<<"original kinematics"<<endl;
                    cout<<"k:  "<<k[0]<<" "<<k[1]<<" "<<k[2]<<endl;
                    cout<<"l:  "<<l[0]<<" "<<l[1]<<" "<<l[2]<<endl;
                    cout<<"q:  "<<q[0]<<" "<<q[1]<<" "<<q[2]<<endl;
                    cout<<"kinematics"<<endl;
                    cout<<"k:  "<<ks[0]<<" "<<ks[1]<<" "<<ks[2]<<" "<<ks[3]<<endl;
                    cout<<"l:  "<<ls[0]<<" "<<ls[1]<<" "<<ls[2]<<" "<<ls[3]<<endl;
                    cout<<"q:  "<<qs[0]<<" "<<qs[1]<<" "<<qs[2]<<" "<<qs[3]<<endl;
                    cout<<"important momenta"<<endl;
                    cout<<"k:  "<<ks[0]+qs[0]<<" "<<ks[1]+qs[1]<<" "<<ks[2]+qs[2]<<" "<<ks[3]+qs[3]<<endl;
                    cout<<"l:  "<<ls[0]<<" "<<ls[1]<<" "<<ls[2]<<" "<<ls[3]<<endl;
                    cout<<"flow check: "<<-ek+eq<<" "<<my_flow.t_val(ni,nf,q)*norm(vector_plus(k,q))<<endl;
                    cout<<"t value: "<<my_flow.t_val(ni,nf,q)<<endl;
                    cout<<"flow jacobian: "<<jac_flow<<endl;
                    cout<<"h function: "<<flow_h<<endl;
                    cout<<"back 2 back jacobian: "<<obs_res.jac<<endl;
                    cout<<"integrand value: "<<integrand_val<<endl;
                    cout<<"product of energies: "<<abs(el*ek*(-ek+eq)*eq)<<endl;
                    cout<<"product of propagators: "<<SP(vector_plus(ks,ls),vector_plus(ks,ls))*SP(vector_plus(ks,vector_plus(ls,qs)),vector_plus(ks,vector_plus(ls,qs)))<<endl;
                    cout<<"final result: "<<jac_flow*flow_h*obs_res.jac*integrand_val<<endl;
                };

                obs_res.jac=obs_res.jac*jac_flow*flow_h;

                dreal colour_average=group_average(obs_res.spin1,obs_res.spin2);

                obs_res.eval=colour_average*integrand_val;

                return obs_res;

            };

            if(comp==3){

                if(DEBUG==1){cout<<"-------------------COMPONENT: "<<comp<<" -------------------------"<<endl;};
                
                vector<vector<dreal>> ni={vector_swap(vector_plus(k,l)),vector_plus(k,vector_plus(l,q))};
                vector<vector<dreal>> nf={};

                if(my_flow.t_val(ni,nf,q)==999){
                    observable res;
                    vector<dreal> j1={0,0,0,0};
                    vector<dreal> j2={0,0,0,0};
                    vector<dreal> pg={0,0,0,0};
                    res.eval=0;
                    res.jac=0;
                    res.j1=j1;
                    res.j2=j2;
                    res.spin1=0;
                    res.spin2=0;
                    res.pg=pg;
                    return res;
                };

                vector<dreal> ks=my_flow.perform_flow(ni,nf,q,k);
                vector<dreal> ls=my_flow.perform_flow(ni,nf,q,l);
                vector<dreal> qs=my_flow.perform_flow(ni,nf,q,q);

                dreal jac_flow=my_flow.jacques(ni,nf,q,3);
                dreal flow_h=my_flow.h(my_flow.t_val(ni,nf,q));

                if(DEBUG==1){
                    cout<<"k:  "<<k[0]<<" "<<k[1]<<" "<<k[2]<<endl;
                    cout<<"l:  "<<l[0]<<" "<<l[1]<<" "<<l[2]<<endl;
                    cout<<"q:  "<<q[0]<<" "<<q[1]<<" "<<q[2]<<endl;
                    cout<<"scaled_momenta"<<endl;
                    cout<<"k:  "<<ks[0]<<" "<<ks[1]<<" "<<ks[2]<<endl;
                    cout<<"l:  "<<ls[0]<<" "<<ls[1]<<" "<<ls[2]<<endl;
                    cout<<"q:  "<<qs[0]<<" "<<qs[1]<<" "<<qs[2]<<endl;
                    cout<<"t_value"<<endl;
                    cout<<my_flow.t_val(ni,nf,q)<<endl;
                    cout<<"scaling jacobian"<<endl;
                    cout<<my_flow.jacques(ni,nf,q,3)<<endl;
                };

                
                dreal eq=sqrt(pow(norm(qs),2.)+pow(m,2));
                dreal el=norm(ls);  
                dreal ekl=norm(vector_plus(ks,ls));

                ks.insert(ks.begin(), el-ekl);
                ls.insert(ls.begin(), -el);
                qs.insert(qs.begin(), eq);


                vector<vector<dreal>> constituents={vector_swap(vector_plus(ks,ls)),vector_plus(ks,vector_plus(ls,qs)),vector_swap(ls)};

                int spins[3]={1,-1,0};

                observable obs_res=my_obs.b2b_sampling_final(constituents, 3, vector_swap(ls), spins, a);

                if(DEBUG==1){
                    vector<dreal> ksp={0.3, 0.1, 0.5};
                    vector<dreal> lsp={0.2, 0.7, 0.1};
                    vector<dreal> qsp={0.1, 0.4, 0.5};

                    dreal ekp=norm(ksp);
                    dreal eqp=sqrt(pow(norm(qsp),2.)+pow(m,2));
                    dreal elp=norm(lsp);  
                    dreal eklp=norm(vector_plus(ksp,lsp));
                    dreal elqp=norm(vector_plus(lsp,qsp));  
                    dreal eklqp=norm(vector_plus(lsp,vector_plus(ksp,qsp)));

                    ksp.insert(ksp.begin(), elp-eklp);
                    lsp.insert(lsp.begin(), -elp);
                    qsp.insert(qsp.begin(), eqp);

                    //NUMERATOR HERE IS OK

                    cout<<"m: "<<m<<endl;
                    cout<<"ksp: "<<ksp[0]<<" "<<ksp[1]<<" "<<ksp[2]<<" "<<ksp[3]<<endl;
                    cout<<"lsp: "<<lsp[0]<<" "<<lsp[1]<<" "<<lsp[2]<<" "<<lsp[3]<<endl;
                    cout<<"qsp: "<<qsp[0]<<" "<<qsp[1]<<" "<<qsp[2]<<" "<<qsp[3]<<endl;
                    cout<<"num: "<<NLO_DT_DrellYan_numerator(vector_swap(vector_plus(ksp,lsp)),vector_swap(ksp),qsp)<<endl;
                    cout<<"den: "<<(16*abs(elp*eklp*(-eklp+eqp)*eqp)*SP(vector_plus(ksp,qsp),vector_plus(ksp,qsp))*SP(ksp,ksp))<<endl;

                    cout<<NLO_DT_DrellYan_numerator(vector_swap(vector_plus(ksp,lsp)),vector_swap(ksp),qsp)/
                    (16*abs(elp*eklp*(-eklp+eqp)*eqp)*SP(vector_plus(ksp,qsp),vector_plus(ksp,qsp))*SP(ksp,ksp))<<endl;
                };


                dreal integrand_val=2*constant*NLO_DT_DrellYan_numerator(vector_swap(vector_plus(ks,ls)),vector_swap(ks),qs)/
                    (16*abs(el*ekl*(-ekl+eq)*eq)*SP(vector_plus(ks,qs),vector_plus(ks,qs))*SP(ks,ks));

                
                if(DEBUG==1){
                    cout<<"############### cut 1 info #############"<<endl;
                    cout<<"original kinematics"<<endl;
                    cout<<"k:  "<<k[0]<<" "<<k[1]<<" "<<k[2]<<endl;
                    cout<<"l:  "<<l[0]<<" "<<l[1]<<" "<<l[2]<<endl;
                    cout<<"q:  "<<q[0]<<" "<<q[1]<<" "<<q[2]<<endl;
                    cout<<"kinematics"<<endl;
                    cout<<"k:  "<<ks[0]<<" "<<ks[1]<<" "<<ks[2]<<" "<<ks[3]<<endl;
                    cout<<"l:  "<<ls[0]<<" "<<ls[1]<<" "<<ls[2]<<" "<<ls[3]<<endl;
                    cout<<"q:  "<<qs[0]<<" "<<qs[1]<<" "<<qs[2]<<" "<<qs[3]<<endl;
                    cout<<"important momenta"<<endl;
                    cout<<"k:  "<<ks[0]+qs[0]<<" "<<ks[1]+qs[1]<<" "<<ks[2]+qs[2]<<" "<<ks[3]+qs[3]<<endl;
                    cout<<"l:  "<<ls[0]<<" "<<ls[1]<<" "<<ls[2]<<" "<<ls[3]<<endl;
                    cout<<"flow check: "<<-ekl+eq<<" "<<my_flow.t_val(ni,nf,q)*norm(vector_plus(k,vector_plus(l,q)))<<endl;
                    cout<<"t value: "<<my_flow.t_val(ni,nf,q)<<endl;
                    cout<<"flow jacobian: "<<jac_flow<<endl;
                    cout<<"h function: "<<flow_h<<endl;
                    cout<<"back 2 back jacobian: "<<obs_res.jac<<endl;
                    cout<<"integrand value: "<<integrand_val<<endl;
                    cout<<"product of energies: "<<abs(el*ekl*(-ekl+eq)*eq)<<endl;
                    cout<<"product of propagators: "<<SP(vector_plus(ks,qs),vector_plus(ks,qs))*SP(ks,ks)<<endl;
                    cout<<"final result: "<<jac_flow*flow_h*obs_res.jac*integrand_val<<endl;
                };

                obs_res.jac=obs_res.jac*jac_flow*flow_h;

                dreal colour_average=group_average(obs_res.spin1,obs_res.spin2);

                obs_res.eval=colour_average*integrand_val;

                return obs_res;

            };

            if(comp==4){

                if(DEBUG==1){cout<<"-------------------COMPONENT: "<<comp<<" -------------------------"<<endl;};
                
                vector<vector<dreal>> ni={vector_swap(l),vector_swap(k),vector_plus(k,vector_plus(l,q))};
                vector<vector<dreal>> nf={};

                if(my_flow.t_val(ni,nf,q)==999){
                    observable res;
                    vector<dreal> j1={0,0,0,0};
                    vector<dreal> j2={0,0,0,0};
                    vector<dreal> pg={0,0,0,0};
                    res.eval=0;
                    res.jac=0;
                    res.j1=j1;
                    res.j2=j2;
                    res.spin1=0;
                    res.spin2=0;
                    res.pg=pg;
                    return res;
                };

                vector<dreal> ks=my_flow.perform_flow(ni,nf,q,k);
                vector<dreal> ls=my_flow.perform_flow(ni,nf,q,l);
                vector<dreal> qs=my_flow.perform_flow(ni,nf,q,q);

                dreal jac_flow=my_flow.jacques(ni,nf,q,3);
                dreal flow_h=my_flow.h(my_flow.t_val(ni,nf,q));

                if(DEBUG==1){
                    cout<<"k:  "<<k[0]<<" "<<k[1]<<" "<<k[2]<<endl;
                    cout<<"l:  "<<l[0]<<" "<<l[1]<<" "<<l[2]<<endl;
                    cout<<"q:  "<<q[0]<<" "<<q[1]<<" "<<q[2]<<endl;
                    cout<<"scaled_momenta"<<endl;
                    cout<<"k:  "<<ks[0]<<" "<<ks[1]<<" "<<ks[2]<<endl;
                    cout<<"l:  "<<ls[0]<<" "<<ls[1]<<" "<<ls[2]<<endl;
                    cout<<"q:  "<<qs[0]<<" "<<qs[1]<<" "<<qs[2]<<endl;
                    cout<<"t_value"<<endl;
                    cout<<my_flow.t_val(ni,nf,q)<<endl;
                    cout<<"scaling jacobian"<<endl;
                    cout<<my_flow.jacques(ni,nf,q,3)<<endl;
                };

                
                dreal eq=sqrt(pow(norm(qs),2.)+pow(m,2));
                dreal el=norm(ls);  
                dreal ek=norm(ks);

                ks.insert(ks.begin(), -ek);
                ls.insert(ls.begin(), -el);
                qs.insert(qs.begin(), eq);

                vector<vector<dreal>> constituents={vector_swap(ls),vector_swap(ks),vector_plus(vector_plus(ks,ls),qs),vector_swap(ls)};

                int spins[4]={0,1,-1,0};

                observable obs_res=my_obs.b2b_sampling_final(constituents, 4, vector_swap(ls), spins, a);

                if(DEBUG==1){
                    vector<dreal> ksp={0.3, 0.1, 0.5};
                    vector<dreal> lsp={0.2, 0.7, 0.1};
                    vector<dreal> qsp={0.1, 0.4, 0.5};

                    dreal ekp=norm(ksp);
                    dreal eqp=sqrt(pow(norm(qsp),2.)+pow(m,2));
                    dreal elp=norm(lsp);  
                    dreal eklp=norm(vector_plus(ksp,lsp));
                    dreal elqp=norm(vector_plus(lsp,qsp));  
                    dreal eklqp=norm(vector_plus(lsp,vector_plus(ksp,qsp)));

                    ksp.insert(ksp.begin(), -ekp);
                    lsp.insert(lsp.begin(), -elp);
                    qsp.insert(qsp.begin(), eqp);

                    //NUMERATOR HERE IS OK

                    cout<<"m: "<<m<<endl;
                    cout<<"ksp: "<<ksp[0]<<" "<<ksp[1]<<" "<<ksp[2]<<" "<<ksp[3]<<endl;
                    cout<<"lsp: "<<lsp[0]<<" "<<lsp[1]<<" "<<lsp[2]<<" "<<lsp[3]<<endl;
                    cout<<"qsp: "<<qsp[0]<<" "<<qsp[1]<<" "<<qsp[2]<<" "<<qsp[3]<<endl;
                    cout<<"num: "<<NLO_DT_DrellYan_numerator(vector_swap(vector_plus(ksp,lsp)),vector_swap(ksp),qsp)<<endl;
                    cout<<"den: "<<(16*abs(elp*ekp*(-ekp-elp+eqp)*eqp)*SP(vector_plus(ksp,lsp),vector_plus(ksp,lsp))*SP(vector_plus(ksp,vector_plus(lsp,qsp)),vector_plus(ksp,vector_plus(lsp,qsp))))<<endl;

                    cout<<NLO_DT_DrellYan_numerator(vector_swap(vector_plus(ksp,lsp)),vector_swap(ksp),qsp)/
                    (16*abs(elp*ekp*(-ekp-elp+eqp)*eqp)*SP(vector_plus(ksp,lsp),vector_plus(ksp,lsp))*SP(vector_plus(ksp,qsp),vector_plus(ksp,qsp)))<<endl;
                };

                dreal integrand_val=2*constant*NLO_DT_DrellYan_numerator(vector_swap(vector_plus(ks,ls)),vector_swap(ks),qs)/
                     (16*abs(el*ek*(-ek-el+eq)*eq)*SP(vector_plus(ks,ls),vector_plus(ks,ls))*SP(vector_plus(ks,qs),vector_plus(ks,qs)));
                
                if(DEBUG==1){
                    cout<<"############### cut 1 info #############"<<endl;
                    cout<<"original kinematics"<<endl;
                    cout<<"k:  "<<k[0]<<" "<<k[1]<<" "<<k[2]<<endl;
                    cout<<"l:  "<<l[0]<<" "<<l[1]<<" "<<l[2]<<endl;
                    cout<<"q:  "<<q[0]<<" "<<q[1]<<" "<<q[2]<<endl;
                    cout<<"kinematics"<<endl;
                    cout<<"k:  "<<ks[0]<<" "<<ks[1]<<" "<<ks[2]<<" "<<ks[3]<<endl;
                    cout<<"l:  "<<ls[0]<<" "<<ls[1]<<" "<<ls[2]<<" "<<ls[3]<<endl;
                    cout<<"q:  "<<qs[0]<<" "<<qs[1]<<" "<<qs[2]<<" "<<qs[3]<<endl;
                    cout<<"important momenta"<<endl;
                    cout<<"k:  "<<ks[0]+qs[0]<<" "<<ks[1]+qs[1]<<" "<<ks[2]+qs[2]<<" "<<ks[3]+qs[3]<<endl;
                    cout<<"l:  "<<ls[0]<<" "<<ls[1]<<" "<<ls[2]<<" "<<ls[3]<<endl;
                    cout<<"flow check: "<<-ek-el+eq<<" "<<my_flow.t_val(ni,nf,q)*norm(vector_plus(k,vector_plus(l,q)))<<endl;
                    cout<<"t value: "<<my_flow.t_val(ni,nf,q)<<endl;
                    cout<<"flow jacobian: "<<jac_flow<<endl;
                    cout<<"h function: "<<flow_h<<endl;
                    cout<<"back 2 back jacobian: "<<obs_res.jac<<endl;
                    cout<<"integrand value: "<<integrand_val<<endl;
                    cout<<"product of energies: "<<abs(el*ek*(-ek-el+eq)*eq)<<endl;
                    cout<<"product of propagators: "<<SP(vector_plus(ks,ls),vector_plus(ks,ls))*SP(vector_plus(ks,qs),vector_plus(ks,qs))<<endl;
                    cout<<"final result: "<<jac_flow*flow_h*obs_res.jac*integrand_val<<endl;
                };

                obs_res.jac=obs_res.jac*jac_flow*flow_h;

                dreal colour_average=group_average(obs_res.spin1,obs_res.spin2);

                obs_res.eval=colour_average*integrand_val;

                return obs_res;

            };

            if(comp==5){

                if(DEBUG==1){cout<<"-------------------COMPONENT: "<<comp<<" -------------------------"<<endl;};
                
                vector<vector<dreal>> ni={vector_swap(vector_plus(k,l)),vector_plus(vector_plus(k,l),q)};
                vector<vector<dreal>> nf={};

                if(my_flow.t_val(ni,nf,q)==999){
                    observable res;
                    vector<dreal> j1={0,0,0,0};
                    vector<dreal> j2={0,0,0,0};
                    vector<dreal> pg={0,0,0,0};
                    res.eval=0;
                    res.jac=0;
                    res.j1=j1;
                    res.j2=j2;
                    res.spin1=0;
                    res.spin2=0;
                    res.pg=pg;
                    return res;
                };

                vector<dreal> ks=my_flow.perform_flow(ni,nf,q,k);
                vector<dreal> ls=my_flow.perform_flow(ni,nf,q,l);
                vector<dreal> qs=my_flow.perform_flow(ni,nf,q,q);

                dreal jac_flow=my_flow.jacques(ni,nf,q,3);
                dreal flow_h=my_flow.h(my_flow.t_val(ni,nf,q));

                if(DEBUG==1){
                    cout<<"k:  "<<k[0]<<" "<<k[1]<<" "<<k[2]<<endl;
                    cout<<"l:  "<<l[0]<<" "<<l[1]<<" "<<l[2]<<endl;
                    cout<<"q:  "<<q[0]<<" "<<q[1]<<" "<<q[2]<<endl;
                    cout<<"scaled_momenta"<<endl;
                    cout<<"k:  "<<ks[0]<<" "<<ks[1]<<" "<<ks[2]<<endl;
                    cout<<"l:  "<<ls[0]<<" "<<ls[1]<<" "<<ls[2]<<endl;
                    cout<<"q:  "<<qs[0]<<" "<<qs[1]<<" "<<qs[2]<<endl;
                    cout<<"t_value"<<endl;
                    cout<<my_flow.t_val(ni,nf,q)<<endl;
                    cout<<"scaling jacobian"<<endl;
                    cout<<my_flow.jacques(ni,nf,q,3)<<endl;
                };

                
                dreal eq=sqrt(pow(norm(qs),2.)+pow(m,2));
                dreal ekl=norm(vector_plus(ks,ls));  
                dreal ek=norm(ks);

                ks.insert(ks.begin(), ek);
                ls.insert(ls.begin(), -ekl-ek);
                qs.insert(qs.begin(), eq);

                vector<vector<dreal>> constituents={vector_swap(vector_plus(ks,ls)),vector_plus(vector_plus(ks,ls),qs),ks,vector_swap(vector_plus(ks,ls))};

                int spins[4]={1,-1,-1,1};

                observable obs_res=my_obs.b2b_sampling_final(constituents, 4, vector_swap(ls), spins, a);

                if(DEBUG==1){
                    vector<dreal> ksp={0.3, 0.1, 0.5};
                    vector<dreal> lsp={0.2, 0.7, 0.1};
                    vector<dreal> qsp={0.1, 0.4, 0.5};

                    dreal ekp=norm(ksp);
                    dreal eqp=sqrt(pow(norm(qsp),2.)+pow(m,2));
                    dreal elp=norm(lsp);  
                    dreal eklp=norm(vector_plus(ksp,lsp));
                    dreal elqp=norm(vector_plus(lsp,qsp));  
                    dreal eklqp=norm(vector_plus(lsp,vector_plus(ksp,qsp)));

                    ksp.insert(ksp.begin(), ekp);
                    lsp.insert(lsp.begin(), -eklp-ekp);
                    qsp.insert(qsp.begin(), eqp);

                    //NUMERATOR HERE IS OK

                    cout<<"m: "<<m<<endl;
                    cout<<"ksp: "<<ksp[0]<<" "<<ksp[1]<<" "<<ksp[2]<<" "<<ksp[3]<<endl;
                    cout<<"lsp: "<<lsp[0]<<" "<<lsp[1]<<" "<<lsp[2]<<" "<<lsp[3]<<endl;
                    cout<<"qsp: "<<qsp[0]<<" "<<qsp[1]<<" "<<qsp[2]<<" "<<qsp[3]<<endl;
                    cout<<"num: "<<NLO_DT_DrellYan_numerator(vector_swap(vector_plus(ksp,lsp)),vector_swap(ksp),qsp)<<endl;
                    cout<<"den: "<<(16*abs(eklp*ekp*(-eklp+eqp)*eqp)*SP(vector_plus(ksp,qsp),vector_plus(ksp,qsp))*SP(lsp,lsp))<<endl;

                    cout<<NLO_DT_DrellYan_numerator(vector_swap(vector_plus(ksp,lsp)),vector_swap(ksp),qsp)/
                    (16*abs(eklp*ekp*(-eklp+eqp)*eqp)*SP(vector_plus(ksp,qsp),vector_plus(ksp,qsp))*SP(lsp,lsp))<<endl;
                };

                dreal integrand_val=2*constant*NLO_DT_DrellYan_numerator(vector_swap(vector_plus(ks,ls)),vector_swap(ks),qs)/
                    (16*abs(ekl*ek*(-ekl+eq)*eq)*SP(vector_plus(ks,qs),vector_plus(ks,qs))*SP(ls,ls));
                
                if(DEBUG==1){
                    cout<<"############### cut 1 info #############"<<endl;
                    cout<<"original kinematics"<<endl;
                    cout<<"k:  "<<k[0]<<" "<<k[1]<<" "<<k[2]<<endl;
                    cout<<"l:  "<<l[0]<<" "<<l[1]<<" "<<l[2]<<endl;
                    cout<<"q:  "<<q[0]<<" "<<q[1]<<" "<<q[2]<<endl;
                    cout<<"kinematics"<<endl;
                    cout<<"k:  "<<ks[0]<<" "<<ks[1]<<" "<<ks[2]<<" "<<ks[3]<<endl;
                    cout<<"l:  "<<ls[0]<<" "<<ls[1]<<" "<<ls[2]<<" "<<ls[3]<<endl;
                    cout<<"q:  "<<qs[0]<<" "<<qs[1]<<" "<<qs[2]<<" "<<qs[3]<<endl;
                    cout<<"important momenta"<<endl;
                    cout<<"k:  "<<ks[0]+qs[0]<<" "<<ks[1]+qs[1]<<" "<<ks[2]+qs[2]<<" "<<ks[3]+qs[3]<<endl;
                    cout<<"l:  "<<ls[0]<<" "<<ls[1]<<" "<<ls[2]<<" "<<ls[3]<<endl;
                    cout<<"flow check: "<<-ekl+eq<<" "<<my_flow.t_val(ni,nf,q)*norm(vector_plus(k,vector_plus(l,q)))<<endl;
                    cout<<"t value: "<<my_flow.t_val(ni,nf,q)<<endl;
                    cout<<"flow jacobian: "<<jac_flow<<endl;
                    cout<<"h function: "<<flow_h<<endl;
                    cout<<"back 2 back jacobian: "<<obs_res.jac<<endl;
                    cout<<"integrand value: "<<integrand_val<<endl;
                    cout<<"product of energies: "<<abs(ekl*ek*(-ekl+eq)*eq)<<endl;
                    cout<<"product of propagators: "<<SP(vector_plus(ks,qs),vector_plus(ks,qs))*SP(ls,ls)<<endl;
                    cout<<"final result: "<<jac_flow*flow_h*obs_res.jac*integrand_val<<endl;
                };

                obs_res.jac=obs_res.jac*jac_flow*flow_h;

                dreal colour_average=group_average(obs_res.spin1,obs_res.spin2);

                obs_res.eval=colour_average*integrand_val;

                return obs_res;

            };

            if(comp==6){

                if(DEBUG==1){cout<<"-------------------COMPONENT: "<<comp<<" -------------------------"<<endl;};
                
                vector<vector<dreal>> ni={vector_swap(k),vector_plus(k,q)};
                vector<vector<dreal>> nf={};

                if(my_flow.t_val(ni,nf,q)==999){
                    observable res;
                    vector<dreal> j1={0,0,0,0};
                    vector<dreal> j2={0,0,0,0};
                    vector<dreal> pg={0,0,0,0};
                    res.eval=0;
                    res.jac=0;
                    res.j1=j1;
                    res.j2=j2;
                    res.spin1=0;
                    res.spin2=0;
                    res.pg=pg;
                    return res;
                };

                vector<dreal> ks=my_flow.perform_flow(ni,nf,q,k);
                vector<dreal> ls=my_flow.perform_flow(ni,nf,q,l);
                vector<dreal> qs=my_flow.perform_flow(ni,nf,q,q);

                dreal jac_flow=my_flow.jacques(ni,nf,q,3);
                dreal flow_h=my_flow.h(my_flow.t_val(ni,nf,q));

                if(DEBUG==1){
                    cout<<"k:  "<<k[0]<<" "<<k[1]<<" "<<k[2]<<endl;
                    cout<<"l:  "<<l[0]<<" "<<l[1]<<" "<<l[2]<<endl;
                    cout<<"q:  "<<q[0]<<" "<<q[1]<<" "<<q[2]<<endl;
                    cout<<"scaled_momenta"<<endl;
                    cout<<"k:  "<<ks[0]<<" "<<ks[1]<<" "<<ks[2]<<endl;
                    cout<<"l:  "<<ls[0]<<" "<<ls[1]<<" "<<ls[2]<<endl;
                    cout<<"q:  "<<qs[0]<<" "<<qs[1]<<" "<<qs[2]<<endl;
                    cout<<"t_value"<<endl;
                    cout<<my_flow.t_val(ni,nf,q)<<endl;
                    cout<<"scaling jacobian"<<endl;
                    cout<<my_flow.jacques(ni,nf,q,3)<<endl;
                };

                
                dreal eq=sqrt(pow(norm(qs),2.)+pow(m,2));
                dreal eklq=norm(vector_plus(vector_plus(ks,ls),qs));  
                dreal ek=norm(ks);

                ks.insert(ks.begin(), -ek);
                ls.insert(ls.begin(), +ek-eq-eklq);
                qs.insert(qs.begin(), eq);

                vector<vector<dreal>> constituents={vector_plus(ks,qs),vector_swap(ks),vector_swap(vector_plus(vector_plus(ks,ls),qs)),vector_plus(ks,qs)};

                int spins[4]={-1,1,1,-1};

                observable obs_res=my_obs.b2b_sampling_final(constituents, 4, vector_swap(ls), spins, a);

                if(DEBUG==1){
                    vector<dreal> ksp={0.3, 0.1, 0.5};
                    vector<dreal> lsp={0.2, 0.7, 0.1};
                    vector<dreal> qsp={0.1, 0.4, 0.5};

                    dreal ekp=norm(ksp);
                    dreal eqp=sqrt(pow(norm(qsp),2.)+pow(m,2));
                    dreal elp=norm(lsp);  
                    dreal eklp=norm(vector_plus(ksp,lsp));
                    dreal elqp=norm(vector_plus(lsp,qsp));  
                    dreal eklqp=norm(vector_plus(lsp,vector_plus(ksp,qsp)));

                    ksp.insert(ksp.begin(), -ekp);
                    lsp.insert(lsp.begin(), ekp-eqp-eklqp);
                    qsp.insert(qsp.begin(), eqp);

                    //NUMERATOR HERE IS OK

                    cout<<"m: "<<m<<endl;
                    cout<<"ksp: "<<ksp[0]<<" "<<ksp[1]<<" "<<ksp[2]<<" "<<ksp[3]<<endl;
                    cout<<"lsp: "<<lsp[0]<<" "<<lsp[1]<<" "<<lsp[2]<<" "<<lsp[3]<<endl;
                    cout<<"qsp: "<<qsp[0]<<" "<<qsp[1]<<" "<<qsp[2]<<" "<<qsp[3]<<endl;
                    cout<<"num: "<<NLO_DT_DrellYan_numerator(vector_swap(vector_plus(ksp,lsp)),vector_swap(ksp),qsp)<<endl;
                    cout<<"den: "<<(16*abs(ekp*(-ekp+eqp)*eklqp*eqp)*SP(vector_plus(ksp,lsp),vector_plus(ksp,lsp))*SP(lsp,lsp))<<endl;

                    cout<<NLO_DT_DrellYan_numerator(vector_swap(vector_plus(ksp,lsp)),vector_swap(ksp),qsp)/
                    (16*abs(ekp*(-ekp+eqp)*eklqp*eqp)*SP(vector_plus(ksp,lsp),vector_plus(ksp,lsp))*SP(lsp,lsp))<<endl;
                };

                dreal integrand_val=2*constant*NLO_DT_DrellYan_numerator(vector_swap(vector_plus(ks,ls)),vector_swap(ks),qs)/
                    (16*abs(ek*(-ek+eq)*eklq*eq)*SP(vector_plus(ks,ls),vector_plus(ks,ls))*SP(ls,ls));
                
                if(DEBUG==1){
                    cout<<"############### cut 1 info #############"<<endl;
                    cout<<"original kinematics"<<endl;
                    cout<<"k:  "<<k[0]<<" "<<k[1]<<" "<<k[2]<<endl;
                    cout<<"l:  "<<l[0]<<" "<<l[1]<<" "<<l[2]<<endl;
                    cout<<"q:  "<<q[0]<<" "<<q[1]<<" "<<q[2]<<endl;
                    cout<<"kinematics"<<endl;
                    cout<<"k:  "<<ks[0]<<" "<<ks[1]<<" "<<ks[2]<<" "<<ks[3]<<endl;
                    cout<<"l:  "<<ls[0]<<" "<<ls[1]<<" "<<ls[2]<<" "<<ls[3]<<endl;
                    cout<<"q:  "<<qs[0]<<" "<<qs[1]<<" "<<qs[2]<<" "<<qs[3]<<endl;
                    cout<<"important momenta"<<endl;
                    cout<<"k:  "<<ks[0]+qs[0]<<" "<<ks[1]+qs[1]<<" "<<ks[2]+qs[2]<<" "<<ks[3]+qs[3]<<endl;
                    cout<<"l:  "<<ls[0]<<" "<<ls[1]<<" "<<ls[2]<<" "<<ls[3]<<endl;
                    cout<<"flow check: "<<-ek+eq<<" "<<my_flow.t_val(ni,nf,q)*norm(vector_plus(k,q))<<endl;
                    cout<<"t value: "<<my_flow.t_val(ni,nf,q)<<endl;
                    cout<<"flow jacobian: "<<jac_flow<<endl;
                    cout<<"h function: "<<flow_h<<endl;
                    cout<<"back 2 back jacobian: "<<obs_res.jac<<endl;
                    cout<<"integrand value: "<<integrand_val<<endl;
                    cout<<"product of energies: "<<abs(ek*(-ek+eq)*eklq*eq)<<endl;
                    cout<<"product of propagators: "<<SP(vector_plus(ks,ls),vector_plus(ks,ls))*SP(ls,ls)<<endl;
                    cout<<"final result: "<<jac_flow*flow_h*obs_res.jac*integrand_val<<endl;
                };

                obs_res.jac=obs_res.jac*jac_flow*flow_h;

                dreal colour_average=group_average(obs_res.spin1,obs_res.spin2);

                obs_res.eval=colour_average*integrand_val;

                return obs_res;

            };

            observable obs_res_f;
            obs_res_f.eval=0;
            obs_res_f.jac=0;
            obs_res_f.j1={0,0,0,0};
            obs_res_f.j2={0,0,0,0};
            obs_res_f.spin1=0;
            obs_res_f.spin2=0;
            obs_res_f.pg={0,0,0,0};


            return obs_res_f;
            
        };

        observable NLO_s_channel_DrellYan_integrand(vector<dreal> k, vector<dreal> l, vector<dreal> q, dreal a, int comp){
            
            if(comp==1){
                
                vector<vector<dreal>> ni={l,vector_plus(vector_minus(k,l),q)};
                vector<vector<dreal>> nf={k};

                if(my_flow.t_val(ni,nf,q)==999){
                    observable res;
                    vector<dreal> j1={0,0,0,0};
                    vector<dreal> j2={0,0,0,0};
                    vector<dreal> pg={0,0,0,0};
                    res.eval=0;
                    res.jac=0;
                    res.j1=j1;
                    res.j2=j2;
                    res.spin1=0;
                    res.spin2=0;
                    res.pg=pg;
                    return res;
                };

                vector<dreal> ks=my_flow.perform_flow(ni,nf,q,k);
                vector<dreal> ls=my_flow.perform_flow(ni,nf,q,l);
                vector<dreal> qs=my_flow.perform_flow(ni,nf,q,q);

                dreal jac_flow=my_flow.jacques(ni,nf,q,3);
                dreal flow_h=my_flow.h(my_flow.t_val(ni,nf,q));

                if(DEBUG==1){
                    cout<<"k:  "<<k[0]<<" "<<k[1]<<" "<<k[2]<<endl;
                    cout<<"l:  "<<l[0]<<" "<<l[1]<<" "<<l[2]<<endl;
                    cout<<"q:  "<<q[0]<<" "<<q[1]<<" "<<q[2]<<endl;
                    cout<<"scaled_momenta"<<endl;
                    cout<<"k:  "<<ks[0]<<" "<<ks[1]<<" "<<ks[2]<<endl;
                    cout<<"l:  "<<ls[0]<<" "<<ls[1]<<" "<<ls[2]<<endl;
                    cout<<"q:  "<<qs[0]<<" "<<qs[1]<<" "<<qs[2]<<endl;
                    cout<<"t_value"<<endl;
                    cout<<my_flow.t_val(ni,nf,q)<<endl;
                    cout<<"scaling jacobian"<<endl;
                    cout<<my_flow.jacques(ni,nf,q,3)<<endl;
                };

                dreal ek=norm(ks);
                dreal eq=sqrt(pow(norm(qs),2.)+pow(m,2));
                dreal el=norm(ls);  

                ks.insert(ks.begin(), ek);
                ls.insert(ls.begin(), el);
                qs.insert(qs.begin(), eq);

                vector<vector<dreal>> constituents={ls,vector_plus(vector_minus(ks,ls),qs)};

                int spins[2]={0,1};

                observable obs_res=my_obs.b2b_sampling_final(constituents, 2, ks, spins, a);;

                if(DEBUG==1){
                    vector<dreal> ksp={0.3, 0.1, 0.5};
                    vector<dreal> lsp={0.2, 0.7, 0.1};
                    vector<dreal> qsp={0.1, 0.4, 0.5};

                    dreal ekp=norm(ksp);
                    dreal eqp=sqrt(pow(norm(qsp),2.)+pow(1,2));
                    dreal elp=norm(lsp);  
                    dreal ekqp=norm(vector_plus(ksp,qsp));
                    dreal elqp=norm(vector_plus(lsp,qsp));  
                    dreal eklqp=norm(vector_plus(lsp,vector_plus(ksp,qsp)));

                    ksp.insert(ksp.begin(), ekp);
                    lsp.insert(lsp.begin(), elp);
                    qsp.insert(qsp.begin(), eqp);

                    cout<<"numerator: "<<NLO_SE_DrellYan_numerator(vector_plus(ksp,qsp),vector_plus(vector_minus(ksp,lsp),qsp),qsp)<<endl;
                    cout<<NLO_SE_DrellYan_numerator(vector_plus(ksp,qsp),vector_plus(vector_minus(ksp,lsp),qsp),qsp)/
                    (16*abs(ekp*elp*(ekp-elp+eqp)*eqp)*SP(vector_plus(ksp,qsp),vector_plus(ksp,qsp))*SP(vector_plus(ksp,qsp),vector_plus(ksp,qsp)))<<endl;
                };
                
                dreal integrand_val=constant*NLO_SE_DrellYan_numerator(vector_plus(ks,qs),vector_plus(vector_minus(ks,ls),qs),qs)/
                    (16*abs(ek*el*(ek-el+eq)*eq)*SP(vector_plus(ks,qs),vector_plus(ks,qs))*SP(vector_plus(ks,qs),vector_plus(ks,qs)));
                
                if(DEBUG==1){
                    cout<<"############### cut 1 info #############"<<endl;
                    cout<<"original kinematics"<<endl;
                    cout<<"k:  "<<k[0]<<" "<<k[1]<<" "<<k[2]<<endl;
                    cout<<"l:  "<<l[0]<<" "<<l[1]<<" "<<l[2]<<endl;
                    cout<<"q:  "<<q[0]<<" "<<q[1]<<" "<<q[2]<<endl;
                    cout<<"kinematics"<<endl;
                    cout<<"k:  "<<ks[0]<<" "<<ks[1]<<" "<<ks[2]<<" "<<ks[3]<<endl;
                    cout<<"l:  "<<ls[0]<<" "<<ls[1]<<" "<<ls[2]<<" "<<ls[3]<<endl;
                    cout<<"q:  "<<qs[0]<<" "<<qs[1]<<" "<<qs[2]<<" "<<qs[3]<<endl;
                    cout<<"important momenta"<<endl;
                    cout<<"k:  "<<ks[0]+qs[0]<<" "<<ks[1]+qs[1]<<" "<<ks[2]+qs[2]<<" "<<ks[3]+qs[3]<<endl;
                    cout<<"l:  "<<ls[0]<<" "<<ls[1]<<" "<<ls[2]<<" "<<ls[3]<<endl;
                    cout<<"t value: "<<my_flow.t_val(ni,nf,q)<<endl;
                    cout<<"flow jacobian: "<<jac_flow<<endl;
                    cout<<"h function: "<<flow_h<<endl;
                    cout<<"back 2 back jacobian: "<<obs_res.jac<<endl;
                    cout<<"integrand value: "<<integrand_val<<endl;
                    cout<<"product of energies: "<<abs(ek*el*(-ek+el+eq)*eq)<<endl;
                    cout<<"product of propagators: "<<SP(vector_plus(ks,qs),vector_plus(ks,qs))*SP(vector_plus(ls,qs),vector_plus(ls,qs))<<endl;
                    cout<<"final result: "<<jac_flow*flow_h*obs_res.jac*integrand_val<<endl;
                };

                obs_res.jac=obs_res.jac*jac_flow*flow_h;

                dreal colour_average=group_average(obs_res.spin1,obs_res.spin2);

                obs_res.eval=colour_average*integrand_val;

                return obs_res;

            };

            observable obs_res_f;
            obs_res_f.eval=0;
            obs_res_f.jac=0;
            obs_res_f.j1={0,0,0,0};
            obs_res_f.j2={0,0,0,0};
            obs_res_f.spin1=0;
            obs_res_f.spin2=0;
            obs_res_f.pg={0,0,0,0};


            return obs_res_f;
            
        };







        
};

extern "C" {

/*dreal mZ=9.118800e+01;
dreal eCM=13000;*/
dreal res_c=0.0;
dreal res_s=0.0;


dreal mZ=9.118800e+01;
dreal eCM=13000;

/*dreal res_c=0.025;
dreal res_s=0.025;*/
/*dreal res_c=1;//0.000012;
//dreal res_c=1;
dreal res_s=1;*/

dreal MUV=91.188;

dreal sigma=2.;


dreal constant_LO=-pow(0.1303037631031056,2)*0.389379*pow(10,9)*3/(3*3*2*2*4*eCM*eCM);


integrands my_integrand(mZ,eCM,sigma,MUV,res_c,res_s);

observable res;
observable_c res_complex;

void set_kinematics(double nM, double neCM){my_integrand.set_kinematics(nM,neCM);};
void set_res(double nres_c, double nres_s){my_integrand.set_res(nres_c,nres_s);};
void set_sigma(double nsigma){my_integrand.set_sigma(nsigma);};
void set_defo_parameters(double lambda, double Mij, double bc_lambda){my_integrand.set_deformation_hyperparameters(lambda, Mij, bc_lambda);};

double LO_scalar_eval(double px, double py, double pz, double qx, double qy, double qz, double a){
    vector<dreal> p={px,py,pz};
    vector<dreal> q={qx,qy,qz};
    return double(my_integrand.LO_scalar(p,q,a));
    };

double LO_DrellYan_eval(double px, double py, double pz, double qx, double qy, double qz, double a){
    vector<dreal> p={px,py,pz};
    vector<dreal> q={qx,qy,qz};
    return double(my_integrand.LO_DrellYan_integrand(p,q,a));
    };

void NLO_dt_DrellYan_eval(double kx, double ky, double kz, double lx, double ly, double lz, double qx, double qy, double qz, double a, int comp){
    vector<dreal> k={kx,ky,kz};
    vector<dreal> l={lx,ly,lz};
    vector<dreal> q={qx,qy,qz};
    res_complex=my_integrand.NLO_dt_DrellYan_integrand(k,l,q,a,comp);
    };

void NLO_se_DrellYan_eval(double kx, double ky, double kz, double lx, double ly, double lz, double qx, double qy, double qz, double a, int comp){
    vector<dreal> k={kx,ky,kz};
    vector<dreal> l={lx,ly,lz};
    vector<dreal> q={qx,qy,qz};
    res_complex=my_integrand.NLO_se_DrellYan_integrand(k,l,q,a,comp);
    };

void NLO_st_channel_DrellYan_eval(double kx, double ky, double kz, double lx, double ly, double lz, double qx, double qy, double qz, double a, int comp){
    vector<dreal> k={kx,ky,kz};
    vector<dreal> l={lx,ly,lz};
    vector<dreal> q={qx,qy,qz};
    res=my_integrand.NLO_st_channel_DrellYan_integrand(k,l,q,a,comp);
    };

void NLO_t_channel_qq_DrellYan_eval(double kx, double ky, double kz, double lx, double ly, double lz, double qx, double qy, double qz, double a, int comp){
    vector<dreal> k={kx,ky,kz};
    vector<dreal> l={lx,ly,lz};
    vector<dreal> q={qx,qy,qz};
    res=my_integrand.NLO_t_channel_qq_DrellYan_integrand(k,l,q,a,comp);
    };

void NLO_t_channel_qg_DrellYan_eval(double kx, double ky, double kz, double lx, double ly, double lz, double qx, double qy, double qz, double a, int comp){
    vector<dreal> k={kx,ky,kz};
    vector<dreal> l={lx,ly,lz};
    vector<dreal> q={qx,qy,qz};
    res=my_integrand.NLO_t_channel_qg_DrellYan_integrand(k,l,q,a,comp);
    };

void NLO_u_channel_DrellYan_eval(double kx, double ky, double kz, double lx, double ly, double lz, double qx, double qy, double qz, double a, int comp){
    vector<dreal> k={kx,ky,kz};
    vector<dreal> l={lx,ly,lz};
    vector<dreal> q={qx,qy,qz};
    res=my_integrand.NLO_u_channel_DrellYan_integrand(k,l,q,a,comp);
    };

void NLO_s_channel_DrellYan_eval(double kx, double ky, double kz, double lx, double ly, double lz, double qx, double qy, double qz, double a, int comp){
    vector<dreal> k={kx,ky,kz};
    vector<dreal> l={lx,ly,lz};
    vector<dreal> q={qx,qy,qz};
    res=my_integrand.NLO_s_channel_DrellYan_integrand(k,l,q,a,comp);
    };

void set_uv_mass(double uvm){
    my_integrand.set_MUV(uvm);
    };

double get_res_eval(int x, int y){
    if(x==0){return double(res.eval);};
    if(x==1){return double(res.jac);};
    if(x==2){return double(res.j1[y]);};
    if(x==3){return double(res.spin1);};
    if(x==4){return double(res.j2[y]);};
    if(x==5){return double(res.spin2);};
    if(x==6){return double(res.pg[y]);};
    return 0;
    };

double get_res_complex_eval(int x, int y, int phase){
    if(x==0){
        if(phase==0){
            return double(res_complex.eval.real());
        }else{
            return double(res_complex.eval.imag());
        }
    };
    if(x==1){
        if(phase==0){
            return double(res_complex.jac.real());
        }else{
            return double(res_complex.jac.imag());
        }
    };
    if(x==2){return double(res_complex.j1[y]);};
    if(x==3){return double(res_complex.spin1);};
    if(x==4){return double(res_complex.j2[y]);};
    if(x==5){return double(res_complex.spin2);};
    if(x==6){return double(res_complex.pg[y]);};
    return 0;
    };

};


int main(){



vector<dreal> k;
k.push_back(0.1);
k.push_back(0.2);
k.push_back(0.3);

vector<dreal> lx;
lx.push_back(0.37);
lx.push_back(0.11);
lx.push_back(-0.44);

vector<dreal> qx;
qx.push_back(-0.22);
qx.push_back(0.2);
qx.push_back(0.77);

deformation_field my_def(1.,0.69,1.);

cout<<my_def.get_deformation_dt(k,lx,qx,1,1)<<" "<<my_def.get_deformation_dt(k,lx,qx,1,2)<<" "<<my_def.get_deformation_dt(k,lx,qx,1,3)<<" "<<endl;
cout<<my_def.jacques_deformation_dt(k,lx,qx,1)<<endl;
//cout<<my_def.jacques_deformation_dt(k,lx,qx,1)<<" "<<endl;


/*
deformation_field my_def(1.,0.69,1.);

cout<<my_def.get_deformation_dt(k,lx,qx)[0]<<" "<<my_def.get_deformation_dt(k,lx,qx)[1]<<" "<<my_def.get_deformation_dt(k,lx,qx)[2]<<" "<<endl;
cout<<my_def.jacques_deformation_dt(k,lx,qx,1)<<" "<<endl;





vector<dreal> kx;
kx.push_back(1);
kx.push_back(0.2);
kx.push_back(0.3);
kx.push_back(0.13);

boost_flow mybf;
rotation_flow myrf;
causal_flow mycf(10.,1.);

dreal t=mybf.t_val(k,1);

vector<dreal> vec=mybf.perform_flow(kx,1,t);
vector<dreal> vecz=mybf.perform_jacques(kx,0.5,1);


observables my_obs(0.1,0.1,1000,90);



vector<vector<dreal>> constituents;

constituents.push_back(kx);
constituents.push_back(lx);
constituents.push_back(qx);

int spins[3]={1,-1,1};

observable res;
res=my_obs.b2b_sampling_final(constituents, 3, kx, spins, 0.5);



integrands my_integrando(91.18, 1000, 1, 1);

vector<dreal> l;
l.push_back(0.63);
l.push_back(0.11);
l.push_back(-0.35);
vector<dreal> q;
q.push_back(1);
q.push_back(0.);
q.push_back(0.);


res=my_integrando.NLO_s_channel_DrellYan_integrand(k,l,q,0.5,1);

cout<<dreal(res.eval)<<endl;

res=my_integrando.NLO_st_channel_DrellYan_integrand(k,l,q,0.5,1);

cout<<dreal(res.eval)<<endl;

res=my_integrando.NLO_u_channel_DrellYan_integrand(k,l,q,0.5,1);

cout<<dreal(res.eval)<<endl;

res=my_integrando.NLO_t_channel_qq_DrellYan_integrand(k,l,q,0.5,1);

cout<<dreal(res.eval)<<endl;

res=my_integrando.NLO_t_channel_qg_DrellYan_integrand(k,l,q,0.5,1);

cout<<dreal(res.eval)<<endl;*/




/*vector<dreal> p={0.1,0.2,0.3};
vector<dreal> q={0.4,0.5,0.1};


boost_flow bf;
rotation_flow rf;

vector<dreal> p12={sqrt(50),4,3,0};


deformation_field my_field(0.99,0.1,0.99);

vector<dreal> kappa=my_field.get_deformation(p,q,1);
dcomp jac=my_field.jacques_deformation(p,q,1,1);

integrands my_int(9.118800e+01,4,1,1,0.3,0.3);

vector<dreal> k={0.1, 0.2, 0.3};
vector<dreal> l={100, 10, 10};
vector<dreal> qn={ 0., 0., 1.};


//vector<dreal> res=my_int.NLO_s_channel_DrellYan_integrand(k,l,qn,0.,1);


vector<dreal> res2=my_int.NLO_se_jj_integrand(k,l,1,1,0);*/

//vector<dreal> res2=my_int.NLO_t_channel_qg_DrellYan_integrand(k,l,qn,0.,1);


/* vector<dreal> np12=bf.perform_flow(p12,1);

vector<dreal> nnp12=bf.perform_flow(np12,2);

for(int i=0; i<4; i++){
    cout<<np12[i]<<" ";
};
cout<<endl;

for(int i=0; i<4; i++){
    cout<<nnp12[i]<<" ";
};
cout<<endl;
*/
/*vector<dreal> p12j={1,2,3,4};

vector<dreal> resj=bf.perform_jacques(p12j, 0.2, 1);

for(int i=0; i<4; i++){
    cout<<resj[i]<<" ";
};
cout<<endl;

cout<<bf.h(0.3)<<endl;*/

/*cout<<"-----"<<endl;

vector<dreal> p12r={sqrt(50),4,3,0.5};

vector<dreal> np12r=rf.perform_flow(p12r,1);

vector<dreal> nnp12r=rf.perform_flow(np12r,2);

for(int i=0; i<4; i++){
    cout<<np12r[i]<<" ";
};
cout<<endl;

for(int i=0; i<4; i++){
    cout<<nnp12r[i]<<" ";
};
cout<<endl;

vector<dreal> p12jr={1,2,3,4};

vector<dreal> resjr=rf.perform_jacques(p12jr, 0.2, 1);

for(int i=0; i<4; i++){
    cout<<resjr[i]<<" ";
};
cout<<endl;

cout<<rf.h(0.3)<<endl;

causal_flow cf(1,1);

vector<vector<dreal>> vi={{4,2.5,0.5},{-1,-0.5,-0.5}};
vector<vector<dreal>> vf={};
vector<dreal> q={5,3,1};

cout<<cf.t_val(vi,vf,q)<<endl;
cout<<cf.jacques(vi,vf,q)<<endl;
cout<<cf.h(0.3)<<endl;*/


/*observables my_obs(0.0, 0.0, 1, 1);

vector<vector<dreal>> constituents={{sqrt(0.04+0.16+0.01),0.2,0.4,0.1},{sqrt(0.25+0.09+0.0),0.5,-0.3,0.0}};

int spins[2]={1,-1};

vector<dreal> my_obs_res=my_obs.b2b_sampling(constituents,2,spins);

cout<<my_obs_res[0]<<" "<<my_obs_res[1]<<" "<<my_obs_res[2]<<" "<<my_obs_res[3]<<" "<<my_obs_res[4]<<" "<<endl;*/

/*
vector<dreal> pn={0.22,0.47,0.111};
vector<dreal> qn={0.65,-0.21,0.38};

vector<dreal> knn={0.2, -0.1, 0.2};
vector<dreal> lnn={0.79, 0.22, 0.413};
vector<dreal> qnn ={0.5, 0.45, 0.37};

integrands my_integrand(1,1,2,1,0.5,0.5);
vector<dreal> res1=my_integrand.NLO_u_channel_DrellYan_integrand(knn,lnn,qnn,0.2,3);


//vector<dreal> res2=my_integrand.NLO_u_channel_DrellYan_integrand(knn,lnn,qnn,0.2,3);

cout<<"res: "<<res1[0]<<" "<<res1[1]<<" "<<res1[2]<<" "<<res1[3]<<" "<<res1[4]<<endl;*/
/*cout<<"rot res: "<<res2[0]<<" "<<res2[1]<<" "<<res2[2]<<" "<<res2[3]<<" "<<res2[4]<<endl;


cout<<"here"<<endl;

vector<vector<dreal>> my_vecs={{1,0,0,1},{1,1/2.,1/2.,1/Sqrt(2)}};

vector<vector<dreal>> rot_vecs = rotator_4d(my_vecs);

cout<<rot_vecs[0][0]<<" "<<rot_vecs[0][1]<<" "<<rot_vecs[0][2]<<" "<<rot_vecs[0][3]<<endl;
cout<<rot_vecs[1][0]<<" "<<rot_vecs[1][1]<<" "<<rot_vecs[1][2]<<" "<<rot_vecs[1][3]<<endl;

vector<dreal> spat_vec1={rot_vecs[0][1],rot_vecs[0][2],rot_vecs[0][3]};
vector<dreal> spat_vec2={rot_vecs[1][1],rot_vecs[1][2],rot_vecs[1][3]};

cout<<norm(spat_vec1)<<" "<<norm(spat_vec2)<<endl;*/

/*vector<dreal> ks={2, 0.3, 0.1, 0.5};
vector<dreal> ls={1, 0.3, 0.7, 0.1};
vector<dreal> qs={0.9, 0.3, 0.1, 0.5};

cout<<my_integrand.NLO_DT_DrellYan_numerator(ks,vector_swap(vector_plus(ls,qs)),vector_swap(qs))<<endl;
vector<dreal> mys=my_integrand.NLO_st_channel_DrellYan_integrand(ks,vector_swap(vector_plus(ls,qs)),vector_swap(qs),1);*/

};