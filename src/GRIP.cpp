#include "GRIP.h"

//using namespace Mongoose;


CompRow_Mat_double ColumnSparsifying(CompRow_Mat_double &inMtx){

    // Parameter for dropping threshold, if <=0 do not sparsify
    if(ColSparsing_thresh <= 0)
        return inMtx;

    int m,i,j,k, *Ap, *Aj ;
    double *Ax;

    m = inMtx.dim(0);
    Ap = inMtx.rowptr_ptr();
    Aj = inMtx.colind_ptr();
    Ax = inMtx.val_ptr();

    int nnz_A = Ap[m];
    int droppedCol=0;

    int *Cp = new int [m + 1];
    int *Cj = new int[nnz_A];
    double *Cx = new double[nnz_A];
    CompRow_Mat_double Amtx;

    for(i=0;i<m+1;i++){
        Cp[i] = 0;
    }

    //Column Sparsing
    if(ColSparsing_thresh > 0){
        double sq = sqrt(m) * ColSparsing_thresh;
        // make it column oriented, to reach columns efficiently
        Amtx = csr_transpose(inMtx); 
        Ap = Amtx.rowptr_ptr();
        Aj = Amtx.colind_ptr();
        Ax = Amtx.val_ptr();

        for(i=0;i<m;i++){
            if( (Ap[i+1] - Ap[i] ) > sq)
            {        
               droppedCol++;
               int todroppedNNZ = (Ap[i+1] - Ap[i]) - sq;
               int nnz_col = Ap[i+1] - Ap[i]; 
               
               double *abs_tmp = new double[nnz_col];
               for(int qq = 0 ; qq < nnz_col; qq++){
                  abs_tmp[qq] = abs(Ax[Ap[i]+qq]);
               }               
               
               priority_queue <double > pq( abs_tmp, abs_tmp+ nnz_col);

               int jj;
               for(jj=0; jj < nnz_col-todroppedNNZ;jj++){
                  pq.pop();
               }

               //cout << "keep values greater than " <<pq.top() << " " << jj << " todroppedNNZ " << todroppedNNZ << " total " << nnz_col << endl;            
               int kept = nnz_col-todroppedNNZ;
               Cp[i+1] = Cp[i];
               for(j=Ap[i];j<Ap[i+1];j++){
                    if((abs(Ax[j]) >= pq.top()) && (kept > 0))
                    {
                        k = Cp[i+1]++;
                        Cj[k] = Aj[j];
                        Cx[k] = Ax[j];
                        kept--;
                    }
                    else{
                       todroppedNNZ--;
                    }
                }
                delete abs_tmp;          
             }                 
             else{         
            // for already sparse columns, copy original column without sparsifying
            // ToDo: this can be improved
                Cp[i+1] = Cp[i];
                for(j=Ap[i];j<Ap[i+1];j++){
                    k = Cp[i+1]++;
                    Cj[k] = Aj[j];
                    Cx[k] = Ax[j];
                }
            }
        }
    }

    int myid;
    MPI_Comm_rank(MPI_COMM_WORLD,&myid);
    if(myid==0) cout << "After columns sparsified new nnz: " << Cp[m] << " DroppedCol: "<< droppedCol << " " <<  "ColSparsing factor:"<< ColSparsing_thresh  << endl;

   //Construct the matrix
    CompRow_Mat_double E (m , m, Cp[m] , Cx, Cp, Cj);

    delete Cp, Cj,Cx;  

   //transpose to convert orientation column to row    
    return  csr_transpose(E); //CSC to CSR    
}


void GRIP_Partitioning(CompRow_Mat_double &Ddrop,idx_t nParts, idx_t* part, double imb)
{
        CompRow_Mat_double AAT;
        CompRow_Mat_double AT = csr_transpose(Ddrop);

        double tt2 = MPI_Wtime();
        AAT = spmm_overlap(Ddrop,AT);
        double tt3 = MPI_Wtime();

        int nnz_AAT = AAT.NumNonzeros();
        cout<< "After SPMM nnz on AA^T " <<  nnz_AAT << " SPMM Time " << tt3-tt2 << endl; ;

        int i,k,j;
        int *jiro = AAT.rowptr_ptr();
        int *jjco = AAT.colind_ptr();
        double *jvalo = AAT.val_ptr();
        int n_o = Ddrop.dim(0);

        idx_t *val_;
        idx_t *adjncy;
        idx_t objval;
        idx_t *ewgt, *vwgt;
        idx_t *xadj  = new idx_t [n_o + 1];

        val_ = new idx_t[nnz_AAT];

         // for enabling edge and vertex weights
         int metis_vwght = 0; // 0: balance on #rows
         int metis_ewght = 1;
                  double totalScaleTime=MPI_Wtime();

         double valMax = -DBL_MAX;
         double valMin = DBL_MAX;	
         for(int i = 0; i < nnz_AAT; i++){
            if(jvalo[i] < 0){
               cout << "ERROR val_i < 0 " << endl; exit(0); }
            if(jvalo[i] != 0){	
               if(jvalo[i]< valMin)
                  valMin=jvalo[i];
               else if(jvalo[i] > valMax)
                  valMax = jvalo[i];
            }
         }	
         double scaleRatio = abs((metis_scaleMax) / (valMax-valMin));
         
	  int metis_epsilon_dropped = 0;
         // for casting floating number to integer
         for(int i = 0; i < nnz_AAT; i++){
            if(jvalo[i]==0) val_[i] = 0;
            else{
         // assuming jvalo is taken its absolute value during SPMM
               double newVal = scaleRatio * abs(( jvalo[i] - valMin ));
               //~ double newVal = scaleRatio *  jvalo[i] ;
               if(newVal < metis_epsilon){
                  val_[i] = 0;
                  metis_epsilon_dropped++;
               }
               else
                  val_[i]  = 1 + floor(newVal);
            }
         }
	 std::cout<<"GRIP scale time : "<<MPI_Wtime()-totalScaleTime<<" sn."<<std::endl; 
         
         //int *Ap = A.rowptr_ptr();
         if(metis_vwght != 0)
         {
            //vwgt = new idx_t[n_o+1];
            for(int i = 0 ; i < n_o; i++){
             //   vwgt[i] = (int)(sqrt(Ap[i+1] - Ap[i])/2)+1;
                if (vwgt[i] == 0 )  cout << " ERROR vertex weight= 0 " << endl;
            }
        }
        else
        {
            vwgt=NULL;
        }

        for(int i = 0; i <= n_o; i++){
            xadj[i] = 0;
        }

        adjncy = new idx_t[nnz_AAT];
        double *tmp_val = new double[nnz_AAT];

        int nz_new = 0;

      // eliminate zero and diagonal values then fill tmp_val
        for (k = 0; k < n_o; k++) {
            xadj[k+1] = xadj[k];
            for (j =jiro[k]; j < jiro[k+1]; j++) {
                //if(val_[j] != 0 && jjco[j] != k)
                if(val_[j] != 0)
                {
                    i = xadj[k+1]++;
                    adjncy[i] = jjco[j];
                    if( metis_ewght != 0)
                    {
                        tmp_val[i] = val_[j];
                    }
                }
            }
        }
        nz_new = xadj[k];

        if(metis_ewght != 0)
        {
            ewgt = new idx_t[nz_new];
            for(int i = 0; i < nz_new; i++){
                ewgt[i] = tmp_val[i];
            }
        }
       else {
            ewgt = NULL;
        }

        delete[] tmp_val;
        delete[] val_;

        printf("Graph scaleFactor %d graphEpsilon %e epsilonDopped %d METIS Preprocessing Time %lf\n",metis_scaleMax, metis_epsilon, metis_epsilon_dropped, MPI_Wtime()-tt3);
        //printf("Graph scaleFactor %d graph Epsilon %e METIS Preprocessing Time %lf\n",metis_scaleMax, metis_epsilon, MPI_Wtime()-tt3);

        idx_t options[METIS_NOPTIONS];
        METIS_SetDefaultOptions(options);

		  //~ int metis_recursive = 0; // for calling METIS_PartGraphRecursive (optional)
//
        //int metis_ufactor = dcntl[Controls::part_imbalance] * 1000; // 0.5*1000= 500 (50% imbalance)

        options[METIS_OPTION_SEED]    = metis_seed;
  //      options[METIS_OPTION_UFACTOR] = imb * 1000; // 0.5*1000= 500 (50% imbalance)

        //options[METIS_OPTION_DBGLVL]  = 17;

        idx_t nVertices = n_o;
        idx_t nWeights = 1;

         cout << "Launching METIS " << " Vwgt: " << metis_vwght <<  " Ewgt: " << metis_ewght  <<" UFAC: " << imb  << " Recursive: " << metis_recursive << " part: " << nParts << " nVertices "<< nVertices << " nEdges " << xadj[nVertices] << endl;
        //cout << "Launching METIS " << " Vwgt: " << metis_vwght <<  " Ewgt: " << metis_ewght  <<" UFAC: " << imb * 1000 << " Recursive: " << metis_recursive << " part: " << nParts << endl;

        double t2 = MPI_Wtime();
        int ret;
 	   
        // for calling METIS_PartGraphRecursive (optional)
        if(metis_recursive == 0){
            ret = METIS_PartGraphKway(&nVertices,&nWeights, xadj, adjncy,  vwgt, NULL, ewgt, &nParts, NULL,NULL, options, &objval, part);
        }else{
            ret = METIS_PartGraphRecursive(&nVertices,&nWeights, xadj, adjncy,  vwgt, NULL, ewgt, &nParts, NULL,NULL, options, &objval, part);
        }
   
        cout << "Done with Numerically Aware Partitioning, time " << setprecision(6) << MPI_Wtime()-t2 << " cut " << objval << " return " << ret << endl;
       cout<<ret<<"-"<<nParts<<"--"<<endl;
        nParts = nParts>>1 ;
        cout<<nParts<<"--"<<endl;
        if(objval < 0){
           cout << "*** WARNING ***" << endl << "--- cutsize is negative, integer overflow problem. reduce metis_scaleMax " << endl;
           cout << " *** " << endl;
        }
        
        delete[] xadj;
        delete[] adjncy;
     }
     
    /*void GRIP_Partitioning2(CompRow_Mat_double &Ddrop,idx_t nParts, idx_t* part, double imb) {
        CompRow_Mat_double AAT;
        CompRow_Mat_double AT = csr_transpose(Ddrop);

        double tt2 = MPI_Wtime();
        AAT = spmm_overlap(Ddrop,AT);
        double tt3 = MPI_Wtime();

        int nnz_AAT = AAT.NumNonzeros();
        cout<< "After SPMM nnz on AA^T " <<  nnz_AAT << " SPMM Time " << tt3-tt2 << endl; ;

        int i,k,j;
        int *jiro = AAT.rowptr_ptr();
        int *jjco = AAT.colind_ptr();
        double *jvalo = AAT.val_ptr();
        int n_o = Ddrop.dim(0);

        double *val_;
        Int *adjncy;
        idx_t objval;
        double *ewgt, *vwgt;
         Int *xadj  = new Int [n_o + 1];

        val_ = new double[nnz_AAT];

         // for enabling edge and vertex weights
         int metis_vwght = 0; // 0: balance on #rows
         int metis_ewght = 1;
         
         double valMax = -DBL_MAX;
         double valMin = DBL_MAX;	
         for(int i = 0; i < nnz_AAT; i++){
            if(jvalo[i] < 0){
               cout << "ERROR val_i < 0 " << endl; exit(0); }
            if(jvalo[i] != 0){	
               if(jvalo[i]< valMin)
                  valMin=jvalo[i];
               else if(jvalo[i] > valMax)
                  valMax = jvalo[i];
            }
         }	
         double scaleRatio = abs((metis_scaleMax) / (valMax-valMin));

         int metis_epsilon_dropped = 0;
         // for casting floating number to integer
         for(int i = 0; i < nnz_AAT; i++){
            if(jvalo[i]==0) val_[i] = 0;
            else{
         // assuming jvalo is taken its absolute value during SPMM
               double newVal = scaleRatio * abs(( jvalo[i] - valMin ));
               //~ double newVal = scaleRatio *  jvalo[i] ;
               if(newVal < metis_epsilon){
                  val_[i] = 0;
                  metis_epsilon_dropped++;
               }
               else
                  val_[i]  = 1 + floor(newVal);
            }
         }
           
         //int *Ap = A.rowptr_ptr();
         if(metis_vwght != 0)
         {
            //vwgt = new idx_t[n_o+1];
            for(int i = 0 ; i < n_o; i++){
             //   vwgt[i] = (int)(sqrt(Ap[i+1] - Ap[i])/2)+1;
                if (vwgt[i] == 0 )  cout << " ERROR vertex weight= 0 " << endl;
            }
        }
        else
        {
            vwgt=NULL;
        }

        for(int i = 0; i <= n_o; i++){
            xadj[i] = 0;
        }

        adjncy = new Int[nnz_AAT];
        double *tmp_val = new double[nnz_AAT];

        int nz_new = 0;

      // eliminate zero and diagonal values then fill tmp_val
        for (k = 0; k < n_o; k++) {
            xadj[k+1] = xadj[k];
            for (j =jiro[k]; j < jiro[k+1]; j++) {
                //if(val_[j] != 0 && jjco[j] != k)
                if(val_[j] != 0)
                {
                    i = xadj[k+1]++;
                    adjncy[i] = jjco[j];
                    if( metis_ewght != 0)
                    {
                        tmp_val[i] = val_[j];
                    }
                }
            }
        }
        nz_new = xadj[k];

        if(metis_ewght != 0)
        {
            ewgt = new double[nz_new];
            for(int i = 0; i < nz_new; i++){
                ewgt[i] = tmp_val[i];
            }
        }
       else {
            ewgt = NULL;
        }

        delete[] tmp_val;
        delete[] val_;

        printf("Graph scaleFactor %d graphEpsilon %e epsilonDopped %d METIS Preprocessing Time %lf\n",metis_scaleMax, metis_epsilon, metis_epsilon_dropped, MPI_Wtime()-tt3);
        //printf("Graph scaleFactor %d graph Epsilon %e METIS Preprocessing Time %lf\n",metis_scaleMax, metis_epsilon, MPI_Wtime()-tt3);

        idx_t options[METIS_NOPTIONS];
        METIS_SetDefaultOptions(options);

		  //~ int metis_recursive = 0; // for calling METIS_PartGraphRecursive (optional)

        //int metis_ufactor = dcntl[Controls::part_imbalance] * 1000; // 0.5*1000= 500 (50% imbalance)

        options[METIS_OPTION_SEED]    = metis_seed;
        options[METIS_OPTION_UFACTOR] = imb * 1000; // 0.5*1000= 500 (50% imbalance)

        //options[METIS_OPTION_DBGLVL]  = 17;

        idx_t nVertices = n_o;
        idx_t nWeights = 1;

         cout << "Launching METIS " << " Vwgt: " << metis_vwght <<  " Ewgt: " << metis_ewght  <<" UFAC: " << imb * 1000  << " Recursive: " << metis_recursive << " part: " << nParts << " nVertices "<< nVertices << " nEdges " << xadj[nVertices] << endl;
        //cout << "Launching METIS " << " Vwgt: " << metis_vwght <<  " Ewgt: " << metis_ewght  <<" UFAC: " << imb * 1000 << " Recursive: " << metis_recursive << " part: " << nParts << endl;

        double t2 = MPI_Wtime();
        int ret;
 	   
        // for calling METIS_PartGraphRecursive (optional)
       /* if(metis_recursive == 0){
            ret = METIS_PartGraphKway(&nVertices,&nWeights, xadj, adjncy,  vwgt, NULL, ewgt, &nParts, NULL,NULL, options, &objval, part);
        }else{
            ret = METIS_PartGraphRecursive(&nVertices,&nWeights, xadj, adjncy,  vwgt, NULL, ewgt, &nParts, NULL,NULL, options, &objval, part);
        }
        EdgeCut_Options *option = EdgeCut_Options::create();
        //option->initial_cut_type = InitialEdgeCut_QP;
        Graph *newCreatedGraph = Graph::create(n_o,nnz_AAT,xadj,adjncy,ewgt,NULL);
        EdgeCut *result = edge_cut(newCreatedGraph, option);
        for (int i = 0; i < newCreatedGraph->n; i++)
    {
           part[i] = result-> partition[i];
    }
        cout << "Done with Numerically Aware Partitioning, time " << setprecision(6) << MPI_Wtime()-t2 << " cut " << objval << " return " << ret << endl;
        if(objval < 0){
           cout << "*** WARNING ***" << endl << "--- cutsize is negative, integer overflow problem. reduce metis_scaleMax " << endl;
           cout << " *** " << endl;
        }
        
        delete[] xadj;
        delete[] adjncy;
     }*/
  

