#include "pTherapy.hh"
#include "Material_Database.hh"
#include "jKiss64.hh"

double Uniform()
{
	double	uni;
	u64 x;
	x=randk();
	uni=(double)x;
	uni/=pow(2.,64.);
	return	uni;
}

double Gaussian(double mean,double rms)
{
	double	data;
	double	U1,U2,Disp,Fx;
  do{
  	U1=2.*Uniform()-1.;
  	U2=2.*Uniform()-1.;
  	Fx=U1*U1+U2*U2;
  }while(Fx>=1.);
  Fx=sqrt((-2.*log(Fx))/Fx);
  Disp=U1*Fx;
  data=mean+Disp*rms;
	return	data;
}

int G4Poisson(double mean)
{
	long	number=0;
	int	 border=16;
	double	position,poissonValue,poissonSum;
	double	value,y,t;
	double	limit=2.e9;
	if(mean<=border)
	{
		position=Uniform();
		poissonValue=exp(-mean);
		poissonSum=poissonValue;
		while(poissonSum<=position) 
		{
			number++;
			poissonValue*=mean/number;
			poissonSum+=poissonValue;
		}
	return	number;
	}
	t=sqrt(-2.*log(Uniform()));
	y=2.*pi*Uniform();
	t*=cos(y);
	value=mean+t*sqrt(mean)+.5;
	if(value<=0.)
		return	0;
	if(value>=limit)
		return	(long)limit;
	return	(long)value;
}

struct Vector Gen4Pi()
{
	double	cosTheta,sinTheta,Phi;
	struct Vector Start_dir;
	cosTheta=2.*Uniform()-1.;
	sinTheta=sqrt(1.-pow(cosTheta,2.));
	Phi=2.*pi*Uniform();
	Start_dir.X=sinTheta*cos(Phi);
	Start_dir.Y=sinTheta*sin(Phi);
	Start_dir.Z=cosTheta;
	return	Start_dir;
}

double Min(double Val1,double Val2)
{
	if(Val1<Val2)
		return	(Val1);
	else 		
		return	(Val2);
}

double Max(double Val1,double Val2)
{
	if(Val1>Val2)	
		return	(Val1);
	else 		
		return	(Val2);
}

struct Vector Troncature(struct Vector Vecteur)
{
	Vecteur.X=epsilon*round(Vecteur.X/epsilon);
	Vecteur.Y=epsilon*round(Vecteur.Y/epsilon);
	Vecteur.Z=epsilon*round(Vecteur.Z/epsilon);
	return Vecteur;
}

struct Vector FirstContact()
{
	double	fDx,fDz,fDx1,fDx2,fDz1,fDz2,angle;
	struct	Vector	FContact;
	
	if(bool.source_rotation==false)
		return	Source.Origine;
	
	angle=Source.Angle*deg;
	FContact.Y=Source.Origine.Y;
	
	if(Source.Origine.X<=Target.Dimension.X&&Source.Origine.X>=-Target.Dimension.X)
	{
		if(angle>90.&&angle<270.)
		{	
			fDz=(Source.Origine.Z-Target.Dimension.Z)/-Source.Direction.Z;
			FContact.Z=Target.Dimension.Z;
		}
		else
		{
			fDz=(Source.Origine.Z+Target.Dimension.Z)/-Source.Direction.Z;
			FContact.Z=-Target.Dimension.Z;
		}
		FContact.X=Source.Origine.X+fDz*Source.Direction.X;
		return	FContact;
	}
	if(Source.Origine.Z<=Target.Dimension.Z&&Source.Origine.Z>=-Target.Dimension.Z)
	{
		if(angle>0.&&angle<180.)
		{	
			fDx=(Source.Origine.X+Target.Dimension.X)/-Source.Direction.X;
			FContact.X=-Target.Dimension.X;
		}
		else
		{
			fDx=(Source.Origine.X-Target.Dimension.X)/-Source.Direction.X;
			FContact.X=Target.Dimension.X;
		}
		FContact.Z=Source.Origine.Z+fDx*Source.Direction.Z;
		return	FContact;
	}
	fDx1=(Source.Origine.X+Target.Dimension.X)/-Source.Direction.X;
	fDx2=(Source.Origine.X-Target.Dimension.X)/-Source.Direction.X;
	fDz1=(Source.Origine.Z+Target.Dimension.Z)/-Source.Direction.Z;
	fDz2=(Source.Origine.Z-Target.Dimension.Z)/-Source.Direction.Z;
	fDx=Min(fDx1,fDx2);
	fDz=Min(fDz1,fDz2);
	if(fDx<fDz)
	{
		if(fDz1<fDz2)
			FContact.Z=-Target.Dimension.Z;
		else
			FContact.Z=Target.Dimension.Z;
		FContact.X=Source.Origine.X+fDz*Source.Direction.X;	
	}
	else
	{
		if(fDx1<fDx2)
			FContact.X=-Target.Dimension.X;
		else
			FContact.X=Target.Dimension.X;
		FContact.Z=Source.Origine.Z+fDx*Source.Direction.Z;	

	}
	return	FContact;
}

void NewSource()
{
	double	U1,U2;
	//Cyrcé
	Source.Dimension.X=
	Source.Dimension.Y=.25*cm;
	do{
			U1=2.*Uniform()-1.;
			U2=2.*Uniform()-1.;
	}while(sqrt(U1*U1+U2*U2)>1.);	
	Source.Origine.X=U1*Source.Dimension.X;
	Source.Origine.X-=1.8*mm;
	Source.Origine.Y=U2*Source.Dimension.Y;
}

void Initialize(int particle,double Energy)
{
	int	i;
	double	U1,U2,tmp;

	id=0;
	for(i=0;i<=idsec;i++)
	{
		Particle[i].nature=
		Particle[i].level=0;
		Particle[i].Ekin=
		Particle[i].mass=
		Particle[i].charge=0.;
		Particle[i].Position.X=
		Particle[i].Position.Y=
		Particle[i].Position.Z=
		Particle[i].Momentum.X=
		Particle[i].Momentum.Y=
		Particle[i].Momentum.Z=0.;
	}
	idsec=0;
	switch(particle)
	{
		case proton:
			Particle[0].mass=m_proton;
			Particle[0].charge=1.;
		break;
		case electron:
			Particle[0].mass=electron_mass_c2;
			Particle[0].charge=-1.;
		break;
		case positron:
			Particle[0].mass=electron_mass_c2;
			Particle[0].charge=1.;
		break;
	}		
	Particle[0].nature=particle;
	Particle[0].Ekin=Energy;

	if(bool.homogenity==true)
	{
		U1=2.*Uniform()-1.;
		U2=2.*Uniform()-1.;
	}
	else
	{
		U1=0.;
		U2=0.;
	}
	Source.Origine.X=U1*Source.Dimension.X;
	Source.Origine.Y=U2*Source.Dimension.Y;
	Source.Origine.Z=Source.Dimension.Z;

	if(bool.source_rotation==true)
	{
		Source.Origine.Z=Source.Dimension.Z*cos(Source.Angle)-Source.Origine.X*sin(Source.Angle);
		Source.Origine.X=Source.Dimension.Z*sin(Source.Angle)+Source.Origine.X*cos(Source.Angle);
		Source.Direction.X=sin(Source.Angle);
		Source.Direction.Y=0.;
		Source.Direction.Z=cos(Source.Angle);
		
	}
	else
	{
		Source.Direction.X=0.;
		Source.Direction.Y=0.;
		Source.Direction.Z=1.;
	}
	if(bool.ray_tracing==true&&Particle[id].level==0)
		fprintf(Ray_tracing_file,"%lf %lf %lf %d\n",Source.Origine.X,Source.Origine.Y,Source.Origine.Z,particle);
	
	if(bool.isotropy==true)
	{
		Source.Direction=Gen4Pi();
		// Source.Direction.X=.02*Source.Origine.X;
		// Source.Direction.Y=.02*Source.Origine.Y;
		// Source.Direction.Z=sqrt(1.-Source.Direction.X*Source.Direction.X-Source.Direction.Y*Source.Direction.Y);
	}
	
	// NewSource();
	Particle[0].Position=FirstContact();
	Particle[0].Momentum=Source.Direction;
	
	for(i=0;i<MAX_PRC;i++)
	{
		Processus[i].CrossSection=0.;
		Processus[i].PreviousCS=0.;
		Processus[i].TrueLength=0.;
	}
}

void Parameters()
{
	int particle;
	particle=Particle[id].nature;
	bool.mapping=false;
	bool.multiplescattering=false;
	bool.perte=false;
	bool.fluctuation=false;
	switch(particle)
	{
		case proton:
			if(bool.hmapping==true)
				bool.mapping=true;
			if(bool.hmultiplescattering==true)
				bool.multiplescattering=true;
			if(bool.hperte==true)
				bool.perte=true;
			if(bool.hfluctuation==true)
				bool.fluctuation=true;							
		break;
		case electron:
		case positron:
			if(bool.emapping==true)
				bool.mapping=true;
			if(bool.emultiplescattering==true)
				bool.multiplescattering=true;
			if(bool.eperte==true||bool.ebrem==true)
				bool.perte=true;
			if(bool.efluctuation==true)
				bool.fluctuation=true;							
		break;
	}	
}
 
void GenerateNewParticle(int particle,double Energy,struct Vector Mom)
{
	switch(particle)
	{
		case electron:
			idsec++;
			Particle[idsec].nature=electron;
			Particle[idsec].level=Particle[id].level+1;
			Particle[idsec].Ekin=Energy;
			Particle[idsec].mass=electron_mass_c2;
			Particle[idsec].charge=-1.;
			Particle[idsec].Position=Particle[id].Position;
			Particle[idsec].Momentum=Mom;
		break;
		case photon:
		break;
	}
	if(bool.verbose==true)
	{
		printf("Nouvelle particule : %d; Energie : %lf; Level : %d\n",idsec,Particle[idsec].Ekin,Particle[idsec].level);
		printf("Position : ");
		PrintVector(Particle[idsec].Position);
		printf("Moment : ");
		PrintVector(Particle[idsec].Momentum);
	}
}

int Inside(struct Vector Place)
{
	Place.X-=Volume.Center.X;	Place.Y-=Volume.Center.Y;	Place.Z-=Volume.Center.Z;
	Place=Troncature(Place);
	if(fabs(Place.X)>Volume.Dimension.X||fabs(Place.Y)>Volume.Dimension.Y||fabs(Place.Z)>Volume.Dimension.Z) 
		return	false;
	else 	
		return	true; 
}

void FirstVoxel(struct Vector Place,u64 Time)
{
	int	i,j,k,ind_x,ind_y,ind_z,idVoxel;
	i=0;	j=0;	k=0;
	while((i*Voxel.Dimension.X)<=(Target.Dimension.X+Place.X)&&i<Voxel.bin.x)
		i++;
	if((i-1)*Voxel.Dimension.X==(Target.Dimension.X+Place.X)&&Particle[id].Momentum.X<0.)
		if(i>0)
			i--;
		else
		{	
			Particle[id].Ekin=0.;
			return;
		}	
	while((j*Voxel.Dimension.Y)<=(Target.Dimension.Y+Place.Y)&&j<Voxel.bin.y)
		j++;
	if((j-1)*Voxel.Dimension.Y==(Target.Dimension.Y+Place.Y)&&Particle[id].Momentum.Y<0.)
		if(j>0)
			j--;
		else
		{	
			Particle[id].Ekin=0.;
			return;
		}	
	while((k*Voxel.Dimension.Z)<=(Target.Dimension.Z+Place.Z)&&k<Voxel.bin.z)
		k++;
	if((k-1)*Voxel.Dimension.Z==(Target.Dimension.Z+Place.Z)&&Particle[id].Momentum.Z<0.)
		if(k>0)
			k--;
		else
		{	
			Particle[id].Ekin=0.;
			return;
		}	
	ind_x=i-1;	ind_y=j-1;	ind_z=k-1;
	if(ind_x<0||ind_y<0||ind_z<0||ind_x>Voxel.bin.x||ind_y>Voxel.bin.y||ind_z>Voxel.bin.z)
	{
		printf("Hors des frontières\n");
		printf("Level : %d\n",Particle[id].level);
		PrintVector(Particle[id].Position);
		printf("Time : %d\n",Time);
		exit(0);
		Particle[id].Ekin=0.;
		return;
	}
	Volume.Center.X=(ind_x+.5)*Voxel.Dimension.X-Target.Dimension.X;
	Volume.Center.Y=(ind_y+.5)*Voxel.Dimension.Y-Target.Dimension.Y;
	Volume.Center.Z=(ind_z+.5)*Voxel.Dimension.Z-Target.Dimension.Z;
	Volume.Center=Troncature(Volume.Center);
	Voxel.indice.x=ind_x;	Voxel.indice.y=ind_y;	Voxel.indice.z=ind_z;
	idVoxel=ind_z*Voxel.bin.x*Voxel.bin.y+ind_y*Voxel.bin.x+ind_x;
	id_mat=Vecteur_Voxel_mat[idVoxel];
	if(bool.verbose==true)
		printf("Voxel %d %d %d -> %d\n",ind_x,ind_y,ind_z,id_mat);
}

int NextVoxel(int face_cube)
{
	int	ind_x,ind_y,ind_z,idVoxel;
	ind_x=Voxel.indice.x;	ind_y=Voxel.indice.y;	ind_z=Voxel.indice.z;
	if(bool.phantom==false)
		return	false;
	if(face_cube==0)
		return	true;
	if(face_cube==2)
		if(ind_x+1<Voxel.bin.x)
			ind_x++;
		else
			return	false;
	if(face_cube==5)
		if(ind_x>=1)
			ind_x--;
		else
			return	false;		
	if(face_cube==3)
		if(ind_y+1<Voxel.bin.y)
			ind_y++;
		else
			return	false;
	if(face_cube==4)
		if(ind_y>=1)
			ind_y--;
		else
			return	false;		
	if(face_cube==1)
		if(ind_z+1<Voxel.bin.z)
			ind_z++;
		else
			return	false;
	if(face_cube==6)
		if(ind_z>=1)
			ind_z--;
		else
			return	false;	
	Volume.Center.X=(ind_x+.5)*Voxel.Dimension.X-Target.Dimension.X;
	Volume.Center.Y=(ind_y+.5)*Voxel.Dimension.Y-Target.Dimension.Y;
	Volume.Center.Z=(ind_z+.5)*Voxel.Dimension.Z-Target.Dimension.Z;
	Volume.Center=Troncature(Volume.Center);
	Voxel.indice.x=ind_x;	Voxel.indice.y=ind_y;	Voxel.indice.z=ind_z;
	idVoxel=ind_z*Voxel.bin.x*Voxel.bin.y+ind_y*Voxel.bin.x+ind_x;
	id_mat=Vecteur_Voxel_mat[idVoxel];
	if(bool.verbose==true)
		printf("Voxel %d %d %d -> %d\n",ind_x,ind_y,ind_z,id_mat);
	// printf("%d (%lf,%lf,%lf) (%lf,%lf,%lf) %lf\n",face_cube,Particle[id].Position.X,Particle[id].Position.Y,Particle[id].Position.Z,
	//  				Volume.Center.X,Volume.Center.Y,Volume.Center.Z,Particle[id].Ekin);
	return	true;
}

int Backscattering(int face_cube)
{
	if(face_cube==1&&Particle[id].Momentum.Z<0.)
		return	true;
	if(face_cube==6&&Particle[id].Momentum.Z>0.)
		return	true;
	if(face_cube==2&&Particle[id].Momentum.X<0.)
		return	true;
	if(face_cube==5&&Particle[id].Momentum.X>0.)
		return	true;
	if(face_cube==3&&Particle[id].Momentum.Y<0.)
		return	true;
	if(face_cube==4&&Particle[id].Momentum.Y>0.)
		return	true;
	return false;
}

double GetSafety(struct Vector Place)
{
	double	Safety=0.,safX,safY,safZ;
	Place.X-=Volume.Center.X;	Place.Y-=Volume.Center.Y;	Place.Z-=Volume.Center.Z;
	Place=Troncature(Place);

	safX=Min(fabs(Volume.Dimension.X-Place.X),fabs(Volume.Dimension.X+Place.X));
	safY=Min(fabs(Volume.Dimension.Y-Place.Y),fabs(Volume.Dimension.Y+Place.Y));
	safZ=Min(fabs(Volume.Dimension.Z-Place.Z),fabs(Volume.Dimension.Z+Place.Z));
	Safety=Min(safX,safY);
	Safety=Min(Safety,safZ);
	if(Safety<0.)
	{
		printf("Safety négatif\n");
		Safety=0.;
	}
	return	Safety;			
}

double GetStepSide(struct Vector Place,double cutLength,int *side,int bool_fragment)
{
	double	exitStep=cutLength,sidX,sidY,sidZ;
	*side=0;
	Place.X-=Volume.Center.X;	Place.Y-=Volume.Center.Y;	Place.Z-=Volume.Center.Z;
	Place=Troncature(Place);
	if(Particle[id].Momentum.X>0.)
	{	
		sidX=Volume.Dimension.X-Place.X;
		sidX/=Particle[id].Momentum.X;
		if(sidX>=0.)	
		{
			exitStep=sidX;
			*side=2;
		}
		else
			exitStep=DBL_MAX;
	}
	if(Particle[id].Momentum.X<0.)
	{	
		sidX=Volume.Dimension.X+Place.X;
		sidX/=-Particle[id].Momentum.X;
		if(sidX>=0.)	
		{	
			exitStep=sidX;
			*side=5;
		}
		else
			exitStep=DBL_MAX;	
	}
	if(Particle[id].Momentum.X==0.)
		exitStep=DBL_MAX;

	if(Particle[id].Momentum.Y>0.)
	{	
		sidY=Volume.Dimension.Y-Place.Y;
		sidY/=Particle[id].Momentum.Y;
		if(exitStep>sidY && sidY>=0.)
		{
			*side=3;
			exitStep=sidY;
		}	
	}	
	if(Particle[id].Momentum.Y<0.)
	{	
		sidY=Volume.Dimension.Y+Place.Y;
		sidY/=-Particle[id].Momentum.Y;
		if(exitStep>sidY && sidY>=0.)
		{
			*side=4;
			exitStep=sidY;
		}	
	}	
	if(Particle[id].Momentum.Z>0.)
	{	
		sidZ=Volume.Dimension.Z-Place.Z;
		sidZ/=Particle[id].Momentum.Z;
		if(exitStep>sidZ && sidZ>=0.)
		{	
			*side=1;
			exitStep=sidZ;	
		}	
	}	
	if(Particle[id].Momentum.Z<0.)
	{	
		sidZ=Volume.Dimension.Z+Place.Z;
		sidZ/=-Particle[id].Momentum.Z;
		if(exitStep>sidZ && sidZ>=0.)
		{	
			*side=6;
			exitStep=sidZ;
		}	
	}	

	if(exitStep<0.)			
	{	
		printf("Fragment<0 erreur : Fragment %5.2E, Face %d, Energie %2.2E, Position (%5.2E,%5.2E,%5.2E), Moment (%5.2E,%5.2E,%5.2E)\n",
						exitStep,*side,Particle[id].Ekin,Particle[id].Position.X,Particle[id].Position.Y,Particle[id].Position.Z,
						Particle[id].Momentum.X,Particle[id].Momentum.Y,Particle[id].Momentum.Z);			
		Particle[id].Ekin=0.*MeV;
		exit(0);
		return;
	}
	
	if(bool_fragment==true&&exitStep>cutLength)
	{	
		printf("Fragment>cut erreur : Cut %5.2E, Fragment %5.2E, Face %d, Energie %2.2E, Position (%5.2E,%5.2E,%5.2E), Moment (%5.2E,%5.2E,%5.2E)\n",
 						cutLength,exitStep,*side,Particle[id].Ekin,Particle[id].Position.X,Particle[id].Position.Y,Particle[id].Position.Z,
						Particle[id].Momentum.X,Particle[id].Momentum.Y,Particle[id].Momentum.Z);
		printf("%d\n",id);
		Particle[id].Ekin=0.*MeV;
		exit(0);
		return;
	}
	return	exitStep;
}

double GetLambda(double Energy,int index)
{
	int	min=0,max=max_ind,mid;
	int	Ekin_inf,Ekin_sup,particle;
	double	Ya,Yb,Xa,Xb;
	double	Valeur;

	while(min<max)
	{
		mid=(min+max)>>1;
		if(Energy>Energy_tab[id_mat][mid])
			min=mid+1;
		else
		max=mid;
	}
	if(min==0)
		min=1;
	if(min==max_ind)
		min=max_ind-1;
	Ekin_inf=min-1;
	Ekin_sup=min;
	Xa=Energy_tab[id_mat][Ekin_inf];
	Xb=Energy_tab[id_mat][Ekin_sup];
	particle=Particle[id].nature;
	switch(index)
	{
		case 1:
			if(particle==electron)
			{
				Ya=eMSC_CS_tab[id_mat][Ekin_inf];
				Yb=eMSC_CS_tab[id_mat][Ekin_sup];
			}
			if(particle==positron)
			{
				Ya=pMSC_CS_tab[id_mat][Ekin_inf];
				Yb=pMSC_CS_tab[id_mat][Ekin_sup];
			}
			if(particle==proton)
			{
				Ya=hMSC_CS_tab[id_mat][Ekin_inf];
				Yb=hMSC_CS_tab[id_mat][Ekin_sup];
			}
		break;
		case 2:
			if(particle==electron)
			{
				Ya=eIoni_CS_tab[id_mat][Ekin_inf];
				Yb=eIoni_CS_tab[id_mat][Ekin_sup];
			}
			if(particle==positron)
			{
				Ya=pIoni_CS_tab[id_mat][Ekin_inf];
				Yb=pIoni_CS_tab[id_mat][Ekin_sup];
			}
			if(particle==proton)
			{
				Ya=hIoni_CS_tab[id_mat][Ekin_inf];
				Yb=hIoni_CS_tab[id_mat][Ekin_sup];
			}
		break;
		case 3:
			if(particle==electron)
			{
				Ya=eBrem_CS_tab[id_mat][Ekin_inf];
				Yb=eBrem_CS_tab[id_mat][Ekin_sup];
			}
			if(particle==positron)
			{
				Ya=pBrem_CS_tab[id_mat][Ekin_inf];
				Yb=pBrem_CS_tab[id_mat][Ekin_sup];
			}	
		break;
		case 4:
			Ya=pAnni_CS_tab[id_mat][Ekin_inf];
			Yb=pAnni_CS_tab[id_mat][Ekin_sup];
		break;
	}
	Valeur=(Energy-Xa)*(Ya-Yb)/(Xa-Xb)+Ya;
	if(Valeur<0.)
		Valeur=0.;
	return	Valeur;
}

double GetRange(double Energy)
{
	int	min=0,max=max_ind,mid;
	int	Ekin_inf,Ekin_sup,particle;
	double	Ya,Yb,Xa,Xb;
	double	Valeur;
	
	while(min<max)
	{
		mid=(min+max)>>1;
		if(Energy>Energy_tab[id_mat][mid])
			min=mid+1;
		else
		max=mid;
	}
	if(min==0)
		min=1;
	if(min==max_ind)
		min=max_ind-1;
	Ekin_inf=min-1;
	Ekin_sup=min;
	Xa=Energy_tab[id_mat][Ekin_inf];
	Xb=Energy_tab[id_mat][Ekin_sup];
	particle=Particle[id].nature;
	if(particle==electron)
	{
		Ya=eRange_tab[id_mat][Ekin_inf];
		Yb=eRange_tab[id_mat][Ekin_sup];
	}
	if(particle==positron)
	{
		Ya=pRange_tab[id_mat][Ekin_inf];
		Yb=pRange_tab[id_mat][Ekin_sup];
	}
	if(particle==proton)
	{
		Ya=hRange_tab[id_mat][Ekin_inf];
		Yb=hRange_tab[id_mat][Ekin_sup];
	}
	Valeur=(Energy-Xa)*(Ya-Yb)/(Xa-Xb)+Ya;
	if(Valeur<0.)
		Valeur=0.;
	return	Valeur;
}

double GetEnergy(double Range)
{
	int	min=0,max=max_ind,mid;
	int	Range_inf,Range_sup,particle;
	double	Ya,Yb,Xa,Xb;
	double	Valeur;

	particle=Particle[id].nature;
	switch(particle)
	{
		case(electron):
			while(min<max)
			{
				mid=(min+max)>>1;
				if(Range>eRange_tab[id_mat][mid])
					min=mid+1;
				else
				max=mid;
			}
			if(min==0)
				min=1;
			if(min==max_ind)
				min=max_ind-1;
			Range_inf=min-1;
			Range_sup=min;
			Xa=eRange_tab[id_mat][Range_inf];
			Xb=eRange_tab[id_mat][Range_sup];
			Ya=eInvRange_tab[id_mat][Range_inf];
			Yb=eInvRange_tab[id_mat][Range_sup];
		break;	
		case(positron):
			while(min<max)
			{
				mid=(min+max)>>1;
				if(Range>pRange_tab[id_mat][mid])
					min=mid+1;
				else
				max=mid;
			}
			if(min==0)
				min=1;
			if(min==max_ind)
				min=max_ind-1;
			Range_inf=min-1;
			Range_sup=min;
			Xa=pRange_tab[id_mat][Range_inf];
			Xb=pRange_tab[id_mat][Range_sup];
			Ya=pInvRange_tab[id_mat][Range_inf];
			Yb=pInvRange_tab[id_mat][Range_sup];
		break;	
		case(proton):
			while(min<max)
			{
				mid=(min+max)>>1;
				if(Range>hRange_tab[id_mat][mid])
					min=mid+1;
				else
				max=mid;
			}
			if(min==0)
				min=1;
			if(min==max_ind)
				min=max_ind-1;
			Range_inf=min-1;
			Range_sup=min;
			Xa=hRange_tab[id_mat][Range_inf];
			Xb=hRange_tab[id_mat][Range_sup];
			Ya=hInvRange_tab[id_mat][Range_inf];
			Yb=hInvRange_tab[id_mat][Range_sup];
		break;	
	}
	Valeur=(Range-Xa)*(Ya-Yb)/(Xa-Xb)+Ya;
	if(Valeur<0.)
		Valeur=0.;
	return	Valeur;
}

double GetDedx(double Energy,int index)
{
	int	min=0,max=max_ind,mid;
	int	Ekin_inf,Ekin_sup,particle;
	double	Ya,Yb,Xa,Xb;
	double	Valeur;

	while(min<max)
	{
		mid=(min+max)>>1;
		if(Energy>Energy_tab[id_mat][mid])
			min=mid+1;
		else
		max=mid;
	}
	if(min==0)
		min=1;
	if(min==max_ind)
		min=max_ind-1;
	Ekin_inf=min-1;
	Ekin_sup=min;
	Xa=Energy_tab[id_mat][Ekin_inf];
	Xb=Energy_tab[id_mat][Ekin_sup];
	particle=Particle[id].nature;
	switch(index)
	{
		case 1:
			if(particle==electron)
			{
				Ya=eIoni_DEDX_tab[id_mat][Ekin_inf];
				Yb=eIoni_DEDX_tab[id_mat][Ekin_sup];
			}
			if(particle==positron)
			{
				Ya=pIoni_DEDX_tab[id_mat][Ekin_inf];
				Yb=pIoni_DEDX_tab[id_mat][Ekin_sup];
			}
			if(particle==proton)
			{
				Ya=hIoni_DEDX_tab[id_mat][Ekin_inf];
				Yb=hIoni_DEDX_tab[id_mat][Ekin_sup];
			}
		break;
		case 2:
			if(particle==electron)
			{
				Ya=eBrem_DEDX_tab[id_mat][Ekin_inf];
				Yb=eBrem_DEDX_tab[id_mat][Ekin_sup];
			}
			if(particle==positron)
			{
				Ya=pBrem_DEDX_tab[id_mat][Ekin_inf];
				Yb=pBrem_DEDX_tab[id_mat][Ekin_sup];
			}	
		break;
	}
	Valeur=(Energy-Xa)*(Ya-Yb)/(Xa-Xb)+Ya;
	if(Valeur<0.)
		Valeur=0.;
	return	Valeur;
}

struct Vector RotateUz(struct Vector newUz,struct Vector vector)
{
	double	u1=newUz.X;
	double	u2=newUz.Y;
	double	u3=newUz.Z;
	double	up=u1*u1+u2*u2;
	double	px,py,pz;
	if(up>0)
	{
		up=sqrt(up);
		px=vector.X,py=vector.Y,pz=vector.Z;
		vector.X=(u1*u3*px-u2*py)/up+u1*pz;
		vector.Y=(u2*u3*px+u1*py)/up+u2*pz;
		vector.Z=-up*px+u3*pz;
	}
	else 
		if(u3<0.)
		{
			vector.X=-vector.X;
			vector.Z=-vector.Z;
		}
	return	vector;
}

struct Vector CorrUnit(struct Vector u, struct Vector v,double uMom, double vMom)
{
	double	r;
	struct	Vector	Final;

	Final.X=u.X*uMom-v.X*vMom;
	Final.Y=u.Y*uMom-v.Y*vMom;
	Final.Z=u.Z*uMom-v.Z*vMom;
	r=Final.X*Final.X+Final.Y*Final.Y+Final.Z*Final.Z;
	if(r>0.)
	{
		r=sqrt(Final.X*Final.X+Final.Y*Final.Y+Final.Z*Final.Z);
		Final.X=Final.X/r;
		Final.Y=Final.Y/r;
		Final.Z=Final.Z/r;
	}
	return	Final;
}

double StepFunction(double Range,double alpha,double rho)
{
	double	StepF;
	if(alpha==0.)
		alpha=.2;
	if(rho==0.)
		rho=1.*mm;	
	if(Range<rho)
		return	Range;
	StepF=alpha*Range+rho*(1.-alpha)*(2.-rho/Range);
	if(StepF<rho)
		StepF=rho;
	return	StepF;
}

double VertexLength(double length,double stepLength)
{
	double	vertexLength;
	if(stepLength>length)
		vertexLength=0.;
	else
		vertexLength=length-stepLength;
	return	vertexLength;
}

void VoxelValidity()
{
	if(fmod(Volume.Dimension.X*2.,Voxel.Dimension.X)>1.E-2)
	{
		Voxel.bin.x=floor(Volume.Dimension.X*2./Voxel.Dimension.X)+1.;
		Volume.Dimension.X=(Voxel.bin.x*Voxel.Dimension.X)/2.;
		printf("Nouvelle dimension sur l'axe x : %5.2lf mm\n",2.*Volume.Dimension.X/mm); 
	}
	else
		Voxel.bin.x=floor(Volume.Dimension.X*2./Voxel.Dimension.X);

	if(fmod(Volume.Dimension.Y*2.,Voxel.Dimension.Y)>1.E-2)
	{
		Voxel.bin.y=floor(Volume.Dimension.Y*2./Voxel.Dimension.Y)+1.;
		Volume.Dimension.Y=(Voxel.bin.y*Voxel.Dimension.Y)/2.;
		printf("Nouvelle dimension sur l'axe y : %5.2lf mm\n",2.*Volume.Dimension.Y/mm); 
	}
	else
		Voxel.bin.y=floor(Volume.Dimension.Y*2./Voxel.Dimension.Y);

	if(fmod(Volume.Dimension.Z*2.,Voxel.Dimension.Z)>1.E-2)
	{
		Voxel.bin.z=floor(Volume.Dimension.Z*2./Voxel.Dimension.Z)+1.;
		Volume.Dimension.Z=(Voxel.bin.z*Voxel.Dimension.Z)/2.;
		printf("Nouvelle dimension sur l'axe z : %5.2lf mm\n",2.*Volume.Dimension.Z/mm); 
	}
	else
		Voxel.bin.z=floor(Volume.Dimension.Z*2./Voxel.Dimension.Z);
		
	if(Voxel.bin.x>MAX_VXL||Voxel.bin.y>MAX_VXL||Voxel.bin.z>MAX_VXL)
		printf("Attention, nombre de voxels trop important\n");	
}

double MapDose(int initialize,double Dose,struct Vector Depot)
{
	int	i,j,k,ind_x,ind_y,ind_z;
	int	nbVoxel,idVoxel;

	Depot.X-=Volume.Center.X;	Depot.Y-=Volume.Center.Y;	Depot.Z-=Volume.Center.Z;
	Depot=Troncature(Depot);
	if(initialize==true)
	{
		EnergyLost=0.;
		nbVoxel=Voxel.bin.x*Voxel.bin.y*Voxel.bin.z;
		Vecteur_Dose_map=(double *)malloc(sizeof(double)*nbVoxel);
		for(i=0;i<nbVoxel;i++)
			Vecteur_Dose_map[i]=0.;
		return;
	}

	if(bool.mapping==false)
		return;

	if(Dose>0.&&id_mat!=0)
	{
		EnergyLost+=Dose;
		Dose=Dose/(Chara_material[id_mat].density*Volume.Dimension.X*Volume.Dimension.Y*Volume.Dimension.Z*8.)*kgramme/J;
		if(bool.phantom==true)
		{
			ind_x=Voxel.indice.x;	ind_y=Voxel.indice.y;	ind_z=Voxel.indice.z;
			idVoxel=ind_z*Voxel.bin.x*Voxel.bin.y+ind_y*Voxel.bin.x+ind_x;
			Vecteur_Dose_map[idVoxel]+=Dose;
			return;
		}
		i=0;	j=0;	k=0;
		while((i*Voxel.Dimension.X)<=(Volume.Dimension.X+Depot.X)&&i<Voxel.bin.x)
			i++;
		while((j*Voxel.Dimension.Y)<=(Volume.Dimension.Y+Depot.Y)&&j<Voxel.bin.y)
			j++;
		while((k*Voxel.Dimension.Z)<=(Volume.Dimension.Z+Depot.Z)&&k<Voxel.bin.z)
			k++;
		ind_x=i-1;	ind_y=j-1;	ind_z=k-1;
		idVoxel=ind_z*Voxel.bin.x*Voxel.bin.y+ind_y*Voxel.bin.x+ind_x;
		Vecteur_Dose_map[idVoxel]+=Dose;
	}
		return;		
}

void RendDose(int initialize,double Dose,struct Vector Depot)
{
	double slice=1000.;
	int	i,bin;
	double depth,binning;
	depth=Depot.X-Volume.Center.X+Volume.Dimension.X;
	binning=Volume.Dimension.X*2./slice;
	depth=Depot.Y-Volume.Center.Y+Volume.Dimension.Y;	
	binning=Volume.Dimension.Y*2./slice;
	depth=Depot.Z-Volume.Center.Z+Volume.Dimension.Z;
	binning=Volume.Dimension.Z*2./slice;
	if(bool.rendement==false)
		return;
		
	if(initialize==true)
	{
		for(i=0;i<MAX_RNG;i++)
			RendementMapping[i]=0;
		return;
	}
	if(Dose>0.)
	{
		bin=(int)ceil(depth/binning);
		RendementMapping[bin]+=Dose;
	}		
}

void RayTracing(double step)
{
	Particle[id].Position.X+=Particle[id].Momentum.X*step;
	Particle[id].Position.Y+=Particle[id].Momentum.Y*step;
	Particle[id].Position.Z+=Particle[id].Momentum.Z*step;
}

int DetectorTrajectory()
{
	double	travel;
	struct	Vector output_point,output_direction;
	output_point.X=Particle[id].Position.X*cos(Source.Angle)-Particle[id].Position.Z*sin(Source.Angle);
	output_point.Y=Particle[id].Position.Y;
	output_point.Z=Particle[id].Position.X*sin(Source.Angle)+Particle[id].Position.Z*cos(Source.Angle);
	output_direction.X=Particle[id].Momentum.X*cos(Source.Angle)-Particle[id].Momentum.Z*sin(Source.Angle);
	output_direction.Y=Particle[id].Momentum.Y;
	output_direction.Z=Particle[id].Momentum.X*sin(Source.Angle)+Particle[id].Momentum.Z*cos(Source.Angle);
	if(output_direction.Z<-epsilon)
		return	false;
	travel=fabs(output_point.Z-Detecteur.Dimension.Z);
	travel/=output_direction.Z;
	output_point.X+=travel*output_direction.X;
	output_point.Y+=travel*output_direction.Y;
	output_point.Z+=travel*output_direction.Z;
	Particle[id].Position.X=output_point.X*cos(Source.Angle)+output_point.Z*sin(Source.Angle);
	Particle[id].Position.Y=output_point.Y;
	Particle[id].Position.Z=-output_point.X*sin(Source.Angle)+output_point.Z*cos(Source.Angle);
	if(bool.ray_tracing==true&&Particle[id].level==0)
		fprintf(Ray_tracing_file,"%lf %lf %lf %d\n",Particle[id].Position.X,Particle[id].Position.Y,Particle[id].Position.Z,Particle[id].nature);
	if(fabs(output_point.X)<Detecteur.Dimension.X&&fabs(output_point.Y)<Detecteur.Dimension.Y)
		return	true;
	else
		return	false;	
}

void AngularCorrection()
{
	struct	Vector tampon;
	tampon.X=Source.Origine.X*cos(Source.Angle)-Source.Origine.Z*sin(Source.Angle);
	tampon.Y=Source.Origine.Y;
	tampon.Z=Source.Origine.X*sin(Source.Angle)+Source.Origine.Z*cos(Source.Angle);
	Source.Origine=tampon;
	tampon.X=Particle[id].Position.X*cos(Source.Angle)-Particle[id].Position.Z*sin(Source.Angle);
	tampon.Y=Particle[id].Position.Y;
	tampon.Z=Particle[id].Position.X*sin(Source.Angle)+Particle[id].Position.Z*cos(Source.Angle);
	Particle[id].Position=tampon;
	tampon.X=Particle[id].Momentum.X*cos(Source.Angle)-Particle[id].Momentum.Z*sin(Source.Angle);
	tampon.Y=Particle[id].Momentum.Y;
	tampon.Z=Particle[id].Momentum.X*sin(Source.Angle)+Particle[id].Momentum.Z*cos(Source.Angle);
	Particle[id].Momentum=tampon;
}

double Theta()
{
	double r,theta;
	
	r=sqrt(Particle[id].Momentum.X*Particle[id].Momentum.X
				+Particle[id].Momentum.Y*Particle[id].Momentum.Y
				+Particle[id].Momentum.Z*Particle[id].Momentum.Z);
	theta=acos(Particle[id].Momentum.Z/r)*deg;
	return	theta;
}

double Phi()
{
	double r,phi;
	
	r=sqrt(Particle[id].Momentum.X*Particle[id].Momentum.X
				+Particle[id].Momentum.Y*Particle[id].Momentum.Y);
	if(r!=0.)
	{
		phi=acos(Particle[id].Momentum.X/r)*deg;	
		if(Particle[id].Momentum.Y<0.)	phi=-phi;
	}
	else
		phi=0.;
	return	phi;	
}

double Tmax()														 
{
	double	tau,ratio,tmax;
	tau=Particle[id].Ekin/Particle[id].mass;
	ratio=electron_mass_c2/Particle[id].mass;
	tmax=2.*electron_mass_c2*tau*(tau+2.)/(1.+2.*(tau+1.)*ratio+ratio*ratio);
	if(tmax>MaxKinEnergy)	
		tmax=MaxKinEnergy;
	return	tmax;
}

double DensCorrection(double x)
{
	double	y=0.;
	
	if(x<Chara_material[id_mat].X0_density)
	{
		if(Chara_material[id_mat].D0_density>0.)
			y=Chara_material[id_mat].D0_density*pow(10.,2.*(x-Chara_material[id_mat].X0_density));
	}
	else
		if(x>=Chara_material[id_mat].X1_density)
			y=2.*log(10.)*x-Chara_material[id_mat].C_density;	
		else
			y=2.*log(10.)*x-Chara_material[id_mat].C_density+Chara_material[id_mat].A_density
			*pow(Chara_material[id_mat].X1_density-x,Chara_material[id_mat].M_density);
	return	y;			
}

double HighOrderCorrections()
{
	int	i,k,Z;
	double	HOC=0.,BarkasTerm=0.,BlochTerm=0.,MottTerm=0.;
	double	tau,gamma,beta,beta2,bg2,ba2,b,W,X,val,y2,term,j,del;

	tau=Particle[id].Ekin/Particle[id].mass;
	gamma=tau+1.;
	bg2=tau*(tau+2.);
	beta2=bg2/(gamma*gamma);
	beta=sqrt(beta2);
	ba2=beta2/(fine_struct*fine_struct);
	//Barkas
	for(i=0;i<Chara_material[id_mat].Nbr_elmt;i++) 
	{
		Z=Chara_material[id_mat].Comp_material[i].Z;
		if(Z==47) 
			BarkasTerm+=Chara_material[id_mat].Comp_material[i].D*Chara_material[id_mat].NbOfAtomsPerVolume*0.006812*pow(beta,-0.9);
		else
		{ 
			if(Z>=64)
				BarkasTerm+=Chara_material[id_mat].Comp_material[i].D*Chara_material[id_mat].NbOfAtomsPerVolume*0.002833*pow(beta,-1.2);
 			else 
 			{    
				X=ba2/Z;
				b=1.3;
				if(Z<=50)	b=1.35;
				if(Z<=25)	b=1.4;
				if(Z<=17)	b=1.4;
				if(Z<=10)	b=1.8;
				if(Z==1)	b=1.8;
				if(Z==2)	b=.6;
				if(Z==18)	b=1.8;
			}
		}	
		W=b/sqrt(X);
		k=0;
		while(BarkasCorr[k][0]<W&&k<46)
			k++;
		if(k==0)
			val=(W-BarkasCorr[k][0])*(BarkasCorr[k][1]-BarkasCorr[k+1][1])
				/(BarkasCorr[k][0]-BarkasCorr[k+1][0])+BarkasCorr[k][1];
		else
			val=(W-BarkasCorr[k][0])*(BarkasCorr[k][1]-BarkasCorr[k-1][1])
				/(BarkasCorr[k][0]-BarkasCorr[k-1][0])+BarkasCorr[k][1];
		if(W>BarkasCorr[46][0]) 
			val*=BarkasCorr[46][0]/W; 
		BarkasTerm+=val*Chara_material[id_mat].Comp_material[i].D*Chara_material[id_mat].NbOfAtomsPerVolume/(sqrt(Z*X)*X);
	}
	BarkasTerm*=1.29*Particle[id].charge/Chara_material[id_mat].NbOfAtomsPerVolume;
	//Bloch
	y2=Particle[id].charge*Particle[id].charge/ba2;
	term=1./(1.+y2);
	j=1.;
	do{
		j++;
		del=1./(j*(j*j+y2));
		term+=del;
	}while(del>.01*term);
	BlochTerm=-y2*term;
	//Mott
	MottTerm=pi*fine_struct*beta*Particle[id].charge;
	//Sum
	HOC=(2.*(BarkasTerm+BlochTerm)+MottTerm);
	return	HOC;
}

double ThetaValue(int KL,int Z)
{
	int	i=0;
	double	Value,Xa,Xb,Ya,Yb,Z_double;
	Z_double=(double)Z;
	switch(KL)
	{
		case 1:
			if(Z_double<ThetaK[0][0])
			{
				Xa=ThetaK[0][0];
				Xb=ThetaK[1][0];
				Ya=ThetaK[0][1];
				Yb=ThetaK[1][1];
			}	
			else
			{
				do{
					i++;
				}while(ThetaK[i][0]<Z_double&&i<35);
				Xa=ThetaK[i-1][0];
				Xb=ThetaK[i][0];
				Ya=ThetaK[i-1][1];
				Yb=ThetaK[i][1];
			}		
		break;
		case 2:
			if(Z_double<ThetaL[0][0])
			{
				Xa=ThetaL[0][0];
				Xb=ThetaL[1][0];
				Ya=ThetaL[0][1];
				Yb=ThetaL[1][1];
			}	
			else
			{
				do{
					i++;
				}while(ThetaL[i][0]<Z_double&&i<36);
				Xa=ThetaL[i-1][0];
				Xb=ThetaL[i][0];
				Ya=ThetaL[i-1][1];
				Yb=ThetaL[i][1];
			}		
		break;

	}
	Value=(Z_double-Xa)*((Ya-Yb)/(Xa-Xb))+Ya;
	return	Value;
}

double Value(double xv,double x1,double x2,double y1, double y2)
{
	return	y1+(y2-y1)*(xv-x1)/(x2-x1);
}

double Value2(double xv,double yv,double x1,double x2,double y1,double y2,double z11,double z21,double z12,double z22)
{
	return	(z11*(x2-xv)*(y2-yv)+z22*(xv-x1)*(yv-y1)
					+.5*(z12*((x2-xv)*(yv-y1)+(xv-x1)*(y2-yv))
					+z21*((xv-x1)*(y2-yv)+(yv-y1)*(x2-xv))))/((x2-x1)*(y2-y1));
}

int Index(double x,int ind_vector,int n)
{
	int indice=n-1;
	switch(ind_vector)
	{
		case 1:
			do{
				indice--;
			}while(indice>0&&x<TheK[indice]);
		break;
		case 2:
			do{
				indice--;
			}while(indice>0&&x<TheL[indice]);
		break;
		case 3:
			do{
				indice--;
			}while(indice>0&&x<Eta[indice]);
		break;
	}
	return	indice;
}

double KShell(double tet,double eta)
{
	double	KShellCorrection=0.;
	int	itet,ieta,nK,nEtaK;
	double	x,y;

	nK=20;
	nEtaK=29;
	x=tet;
	itet=0;
	ieta=0;
	if(tet<TheK[0])
		x=TheK[0];
	else
		if(tet>TheK[nK-1])
		{
			x=TheK[nK-1];
			itet=nK-2;
		}
		else
			itet=Index(x,1,nK);

	if(eta>=Eta[nEtaK-1])
		KShellCorrection=(Value(x,TheK[itet],TheK[itet+1],UK[itet],UK[itet])
										 +Value(x,TheK[itet],TheK[itet+1],VK[itet],VK[itet])/eta
										 +Value(x,TheK[itet],TheK[itet+1],ZK[itet],ZK[itet])/(eta*eta))/eta;
	else
	{
		y=eta;
		if(eta<Eta[0])
			y=Eta[0];
		else
			ieta=Index(y,3,nEtaK);
		KShellCorrection=Value2(x,y,TheK[itet],TheK[itet+1],Eta[ieta],Eta[ieta+1],
										 CK[itet][ieta],CK[itet+1][ieta],CK[itet][ieta+1],CK[itet+1][ieta+1]);
	}
	return	KShellCorrection;
}

double LShell(double tet,double eta)
{
	double	LShellCorrection=0.;
	int	itet,ieta,nL,nEtaL;
	double	x,y;

	nL=26;
	nEtaL=28;
	x=tet;
	itet=0;
	ieta=0;
	if(tet<TheL[0])
		x=TheL[0];
	else
		if(tet>TheL[nL-1])
		{
			x=TheL[nL-1];
			itet=nL-2;
		}
		else
			itet=Index(x,2,nL);

	if(eta>=Eta[nEtaL-1])
		LShellCorrection=(Value(x,TheL[itet],TheL[itet+1],UL[itet],UL[itet])
										 +Value(x,TheL[itet],TheL[itet+1],VL[itet],VL[itet])/eta)/eta;
	else
	{
		y=eta;
		if(eta<Eta[0])
			y=Eta[0];
		else
			ieta=Index(y,3,nEtaL);
		LShellCorrection=Value2(x,y,TheL[itet],TheL[itet+1],Eta[ieta],Eta[ieta+1],
										 CL[itet][ieta],CL[itet+1][ieta],CL[itet][ieta+1],CL[itet+1][ieta+1]);
	}
	return	LShellCorrection;
}

double ShellCorrection()
{
	double	ShellCorrectionTerm=0.;
	int	i,j,Z,ntot,nmax,ne;
	double	res,res0,Z2,Zeff,norm,eshell,tet,eta,f;
	double	alpha2,tau,gamma,bg2,beta2,ba2;

	tau=Particle[id].Ekin/Particle[id].mass;
	gamma=tau+1.;
	bg2=tau*(tau+2.);
	beta2=bg2/(gamma*gamma);
	alpha2=fine_struct*fine_struct;
	ba2=beta2/alpha2;

	for(i=0;i<Chara_material[id_mat].Nbr_elmt;i++)
	{
		res=0.;
		res0=0.;
		Z=Chara_material[id_mat].Comp_material[i].Z;
		Z2=pow((Z-.3),2);
		f=1.;
		if(Z==1)
		{
			f=.5;
			Z2=1.;
		}
		eta=ba2/Z2;
		tet=Z2*(1.+Z2*.25*alpha2);
		if(Z>11)
			tet=ThetaValue(1,Z);
		res0=f*KShell(tet,eta);
		res+=res0;
		if(Z>2)
		{
			Zeff=Z-ZD[10];
			if(Z<10)
				Zeff=Z-ZD[Z];
			Z2=Zeff*Zeff;
			eta=ba2/Z2;
			f=.125;
			tet=ThetaValue(2,Z);
			ntot=Chara_material[id_mat].Comp_material[i].S;
			nmax=(int)Min(4,ntot);
			norm=0.;
			eshell=0.;
			for(j=1;j<nmax;j++)
			{
				ne=Chara_material[id_mat].Comp_material[i].ElecOfShell[j];
				if(Z<=15)
				{
					if(j<3)
						tet=.25*Z2*(1.+5.*Z2*alpha2/16.);
					else
						tet=.25*Z2*(1.+Z2*alpha2/16.);
				}
				norm+=ne;
				eshell+=tet*ne;
				res0=f*ne*LShell(tet,eta);
				res+=res0;
			}
			if(ntot>nmax)
			{
				eshell/=norm;
				if(Z<28)
					res+=f*(Z-10.)*LShell(eshell,HM[Z-11]*eta);
				else
					if(Z<63)
						res+=f*18*LShell(eshell,HM[Z-11]*eta);
					else
						res+=f*18*LShell(eshell,HM[52]*eta);
				if(Z>32)
				{
					if(Z<60)
						res+=f*(Z-28)*LShell(eshell,HN[Z-33]*eta);
					else
						if(Z<63)
							res+=4*LShell(eshell,HN[Z-33]*eta);
						else
							res+=4*LShell(eshell,HN[30]*eta);
					if(Z>60)
						res+=f*(Z-60)*LShell(eshell,150.*eta);
				}	
			}
		}
		ShellCorrectionTerm+=res*Chara_material[id_mat].Comp_material[i].D/Z;
	}
	return	ShellCorrectionTerm;
}			

double hBetheIoniDEDX(double CutEnergyElectron)
{
	int	i,Z;
	double	Dedx=0.;
	double	tau,gamma,bg2,beta2,eexc2,x;
	double	tmax=Tmax(),cutEnergy=Min(CutEnergyElectron,tmax);

	tau=Particle[id].Ekin/Particle[id].mass;
	gamma=tau+1.;
	bg2=tau*(tau+2.);
	beta2=bg2/(gamma*gamma);
	eexc2=Chara_material[id_mat].excitation_energy*Chara_material[id_mat].excitation_energy;
	Dedx=log(2.*electron_mass_c2*bg2*cutEnergy/eexc2)-(1.+cutEnergy/tmax)*beta2;
	Dedx+=pow((.5*cutEnergy/(Particle[id].Ekin+Particle[id].mass)),2);
	x=log(bg2)/(2.*log(10.));
	Dedx-=DensCorrection(x);
	Dedx-=2.*ShellCorrection();
	Dedx+=HighOrderCorrections();
	Dedx*=twopi_mc2_rcl2*Chara_material[id_mat].ElecDensity/beta2;
	if(Dedx<0.)
		Dedx=0.;
	return	Dedx;
}

double eIoniDEDX(double CutEnergyElectron)
{
	double	Dedx=0.;
	double	th=.25*sqrt(Chara_material[id_mat].Zeff)*keV;
	double	lowLimit=.2*keV;
	double	tmax,tkin;
	double	eexc,eexc2,d,x,y;
	double	tau,gamma,gamma2,beta2,bg2;
	double	d2,d3,d4;
	
	tkin=Particle[id].Ekin;
	if(Particle[id].Ekin<th)
		tkin=th;
  tmax=tkin*.5;
 	if(Particle[id].nature==positron)
		tmax=tkin;

  tau=tkin/electron_mass_c2;
  gamma=tau+1.;
  gamma2=gamma*gamma;
 	beta2=1.-1./gamma2;
  bg2=beta2*gamma2;
  eexc=Chara_material[id_mat].excitation_energy/electron_mass_c2;
 	eexc2=eexc*eexc; 
	d=Min(CutEnergyElectron,tmax);
	d/=electron_mass_c2;
	if(Particle[id].nature==electron)
		Dedx=log(2.*(tau+2.)/eexc2)-1.-beta2+log((tau-d)*d)+tau/(tau-d)
				+(.5*d*d+(2.*tau+1.)*log(1.-d/tau))/gamma2;
  else
  {
    d2=d*d*.5;
    d3=d2*d/1.5;
    d4=d3*d*.75;
    y=1./(1.+gamma);
    Dedx=log(2.*(tau+2.)/eexc2)+log(tau*d)
   		  -beta2*(tau+2.*d-y*(3.*d2+y*(d-d3+y*(d2-tau*d3+d4))))/tau;
  } 

  x=log(bg2)/(2.*log(10.));
  Dedx-=DensCorrection(x); 
  Dedx*=twopi_mc2_rcl2*Chara_material[id_mat].ElecDensity/beta2;
  if(Dedx<0.)
  	Dedx=0.;
	if(Particle[id].Ekin<th)
	{
    if (Particle[id].Ekin>=lowLimit) 
    	Dedx*=sqrt(tkin/Particle[id].Ekin);
    else
    	Dedx*=sqrt(tkin*Particle[id].Ekin)/lowLimit;
  }
  return	Dedx;
}

double hElectronicStoppingPower(int Z,double Energy)
{
	double	ElecStopPow=0.;
	int	i;
	double	fac=1.,slow,shigh;
	i=Z-1;
	if(i<0)
		i=0;
	if(i>91)
		i=91;
	Energy=Energy/(m_proton/amu_c2*keV);
	if(Energy<40.&&i==5)
	{
		fac=sqrt(Energy/40.);
		Energy=40.;
	}
	else
		if(Energy<10.)
		{
			fac=sqrt(Energy*.1);
			Energy=10.;
		}
	slow=ESP[i][1]*pow(Energy,.45);
	shigh=log(1.+ESP[i][3]/Energy+ESP[i][4]*Energy)*ESP[i][2]/Energy;
	ElecStopPow=slow*shigh*fac/(slow+shigh);
	if(ElecStopPow<0.)
		ElecStopPow=0.;
	return	ElecStopPow;
}

double hBDEDX(double Energy) 
{
	int	ind_en,n;
	double	eLoss=0.;
	double	xa,xb,ya,yb;
	double	theZieglerFactor=1.E-15*eV*cm2;
	
	if(Energy>2.01*MeV)
	{
		printf("### Warning ###\nOutside the Bragg Model\n");
		return	eLoss;
	}	

	if(pstarE[0]!=0.) 
	{
		if(Energy<pstarT[0]*MeV)
			eLoss=pstarE[0]*sqrt(Energy/(pstarT[0]*MeV));
		else
		{	
			ind_en=0;
			while(pstarT[ind_en]*MeV<Energy&&ind_en<61)
				ind_en++;
			xa=pstarT[ind_en-1];	xb=pstarT[ind_en];	
			ya=pstarE[ind_en-1];	yb=pstarE[ind_en];
			eLoss=(Energy-xa)*(ya-yb)/(xa-xb)+ya;
		}
		eLoss*=Chara_material[id_mat].density*MeV*cm2/gramme;	
	}
	else
	{
		for(n=0;n<Chara_material[id_mat].Nbr_elmt;n++)
			eLoss+=hElectronicStoppingPower(Chara_material[id_mat].Comp_material[n].Z,Energy)
						*Chara_material[id_mat].Comp_material[n].D*Chara_material[id_mat].NbOfAtomsPerVolume
						*theZieglerFactor;
	}	
	return	eLoss;
}

double hBraggDEDX(double cutEnergyElectron)
{
	double	tau,gamma,beta2,bg2,x;
	double	tmax=Tmax(),lowKinEnergy=1.*keV;
	double	energy,Dedx=0.;
	
	energy=Particle[id].Ekin;
	if(energy>lowKinEnergy)
		Dedx=hBDEDX(energy); 
	else
		Dedx=hBDEDX(lowKinEnergy)*sqrt(energy/lowKinEnergy);

	if(tmax>cutEnergyElectron) 
	{
		tau=energy/Particle[id].mass;
		gamma=tau+1.;
		bg2=tau*(tau+2.);
		beta2=bg2/(gamma*gamma);
		x=cutEnergyElectron/tmax;
		Dedx+=(log(x)+(1.-x)*beta2)*twopi_mc2_rcl2*Chara_material[id_mat].ElecDensity/beta2;
	}
	if(Dedx<0.) 
		Dedx=0.;
	Dedx*=Particle[id].charge*Particle[id].charge;

	return	Dedx;
}

double eBremLoss(double Z,double T,double Cut)
{
  int	i,j;
  int	NZ=8,Nloss=11,iz=0;
  double	Loss;
  double	dz,xx,yy,fl,E;
  double	aaa=.414,bbb=.345,ccc=.460,delz=1.e6;
	double	beta=1.0,ksi=2.0,clossh=.254,closslow=1./3.,alosslow=1.;
  double	Tlim=10.*MeV,xlim=1.2;

	for(i=0;i<NZ;i++)
	{
		dz=fabs(Z-ZZ[i]);
		if(dz<delz)
		{
			iz=i;
			delz=dz;
		}
	}
	xx=log10(T);
	fl=1.;
	if(xx<=xlim)
	{
		xx/=xlim;
		yy=1.;
		fl=0.;
		for(j=0;j<Nloss;j++) 
		{
			fl+=yy+coefloss[iz][j];
			yy*=xx;
		}
		if(fl<.00001)
			fl=.00001;
		else 
			if(fl>1.)
				fl=1.;
	}

	E=T+electron_mass_c2;
	Loss=Z*(Z+ksi)*E*E/(T+E)*exp(beta*log(Cut/T))*(2.-clossh*exp(log(Z)/4.));
	if(T<=Tlim)
		Loss/=exp(closslow*log(Tlim/T));
	if(T<=Cut)
		Loss*=exp(alosslow*log(T/Cut));
	Loss*=(aaa+bbb*T/Tlim)/(1.+ccc*T/Tlim);
	Loss*=fl;
	Loss/=N_avogadro;

	return	Loss;
}	

double PositronCorrLossFactor(double Z,double T,double Cut)
{
	double	Factor=0.;
	const double	K=132.9416*eV;
	const double	a1=4.15E-1,a3=2.10E-3,a5=54.E-5;
	double	x,x2,x3,eta,e0;

	x=log(T/(K*Z*Z));
	x2=x*x;
	x3=x2*x;
	eta=.5*atan(a1*x+a3*x3+a5*x3*x2)/pi;
	e0=Cut/T;
	if(e0<1.)
		Factor=exp(log(1.-e0)/eta);
	Factor=eta*(1.-Factor)/e0;
	return	Factor;
}

double eBremDEDX(double cutEnergy)
{
	int	i,n,nn,nmax;
	double	Dedx;
	double	totalEnergy,Z,natom,kp2,kmin,kmax,floss;
	double	vmin,vmax,u,fac,c,v,dv;
	double	thigh=100.*GeV;
	double	cut=Min(cutEnergy,Particle[id].Ekin);
	double	rate,loss;
	double	factorHigh=36./(1450.*GeV);
	double	coef1=-.5;
	double	coef2=2./9.;
	double	lowKinEnergy=0.*eV;
	double	highKinEnergy=1.*GeV;
	double	probsup=1.;
	double	MigdalConstant=elec_radius*hbarc*hbarc*4.*pi/(electron_mass_c2*electron_mass_c2);

	totalEnergy=Particle[id].Ekin+electron_mass_c2;
	Dedx=0.;

	if(Particle[id].Ekin<lowKinEnergy)
		return	0.;

	for(i=0;i<Chara_material[id_mat].Nbr_elmt;i++) 
	{
		Z=Chara_material[id_mat].Comp_material[i].Z;
		natom=Chara_material[id_mat].Comp_material[i].D;
		if(Particle[id].Ekin<=thigh) 
			loss=eBremLoss(Z,Particle[id].Ekin,cut);
		if(Particle[id].nature==positron)
			loss*=PositronCorrLossFactor(Z,Particle[id].Ekin,cut);
		loss*=natom;
		kp2=MigdalConstant*totalEnergy*totalEnergy*Chara_material[id_mat].ElecDensity;

		kmin=1.*eV;
		kmax=cut;
		if(kmax>kmin) 
		{
			floss=0.;
			nmax=100;
			vmin=log(kmin);
			vmax=log(kmax);	
			nn=(int)(nmax*(vmax-vmin)/(log(highKinEnergy)-vmin)) ;
			if(nn>0) 
			{
				dv=(vmax-vmin)/nn;
				v=vmin-dv;
				for(n=0;n<=nn;n++) 
				{
					v+=dv;
					u=exp(v);
					//fac=u*SupressionFunction(material,Particle[id].Ekin,u);	//LPM flag off 
					fac=u*1.;
					fac*=probsup*(u*u/(u*u+kp2))+1.-probsup;
					if((n==0)||(n==nn)) 
						c=.5;
					else
			    	c=1.;
					fac*=c;
					floss+=fac ;
				}
				floss*=dv/(kmax-kmin);
			}
			else
				floss=1.;
			if(floss>1.)
				floss=1.;
			loss*=floss;
		}
	Dedx+=loss;
	}
	if(Dedx<0.)
		Dedx=0.;
	Dedx*=Chara_material[id_mat].NbOfAtomsPerVolume;
	return	Dedx;
}

double LossApproximation(double StepLength)
{
	int i;
	double	range,perteApp;

	range=GetRange(Particle[id].Ekin);
	range-=StepLength;
	perteApp=GetEnergy(range);
	perteApp=Particle[id].Ekin-perteApp;
	return	perteApp;	
}

double hFluctuation(double meanLoss,double CutElectron,double StepLength)
{
	int	i,nb;
	double	LossFluct=0.,minLoss=10.*eV;
	double	minNumberInteractionsBohr,siga,x,tmaxkin,massrate;
	double	tmax,gamma,gamma2,beta2,rate,w,w1,w2,C,fw,lossc,a1,a2,a3,nmaxCont;
	double	sa1,alfa,alfa1,nmean,e1,e2,esmall,emean,p1,p2,p3,sige,sig2e;

	if(meanLoss<minLoss)
		return	meanLoss;
		
	nmaxCont=16.,fw=4.,rate=.55,emean=0.,lossc=0.,sige=0.,sig2e=0.;
	minNumberInteractionsBohr=10.;a1=0.;a2=0.;a3=0.;p1=0.;p2=0.;p3=0.;
	tmax=Min(Tmax(),CutElectron);
	esmall=.5*sqrt(Chara_material[id_mat].E0*Chara_material[id_mat].excitation_energy);
	gamma=Particle[id].Ekin/Particle[id].mass+1.;
	gamma2=gamma*gamma;
	beta2=Particle[id].Ekin/Particle[id].mass*(Particle[id].Ekin/Particle[id].mass+2.)/gamma2;

	if(Particle[id].mass>electron_mass_c2&&meanLoss>=minNumberInteractionsBohr*tmax)
	{
		massrate=electron_mass_c2/Particle[id].mass;
		tmaxkin=2.*electron_mass_c2*beta2*gamma2/(1.+massrate*(2.*gamma+massrate));
		if(tmaxkin<=2.*tmax)   
		{
			siga=(1./beta2-.5)*twopi_mc2_rcl2*tmax*StepLength
					*Chara_material[id_mat].ElecDensity*Particle[id].charge*Particle[id].charge;
			siga=sqrt(siga);
			if(2.*meanLoss<siga) 
				do{
					LossFluct=2.*meanLoss*Uniform();
					x=(LossFluct-meanLoss)/siga;
				}while(1.-.5*x*x<Uniform());
			else 
				do{
					LossFluct=Gaussian(meanLoss,siga);
				}while(LossFluct<0.||LossFluct>2.*meanLoss);
			return	LossFluct;
		}
	}

	if(tmax<=Chara_material[id_mat].E0)
		return	meanLoss;

	if(tmax>Chara_material[id_mat].excitation_energy)
	{
		w2=log(2.*electron_mass_c2*beta2*gamma2)-beta2;
		if(w2>Chara_material[id_mat].excitationLog&&w2>Chara_material[id_mat].E2Log)
		{
			C=meanLoss*(1.-rate)/(w2-Chara_material[id_mat].excitationLog);
			a1=C*Chara_material[id_mat].F1*(w2-Chara_material[id_mat].E1Log)/Chara_material[id_mat].E1;
			a2=C*Chara_material[id_mat].F2*(w2-Chara_material[id_mat].E2Log)/Chara_material[id_mat].E2;
			if(a1<nmaxCont)
			{
				sa1=sqrt(a1);
				if(Uniform()<exp(-sa1))
				{
					e1=esmall;
					a1=meanLoss*(1.-rate)/e1;
					a2=0.;
					e2=Chara_material[id_mat].E2;
				}
				else
				{
					a1=sa1;
					e1=sa1*Chara_material[id_mat].E1;
					e2=Chara_material[id_mat].E2;
				}
			}
			else															
			{
				a1=a1/fw;
				e1=fw*Chara_material[id_mat].E1;
				e2=Chara_material[id_mat].E2;
			}
		}
	}
	w1=tmax/Chara_material[id_mat].E0;
	if(tmax>Chara_material[id_mat].E0)
		a3=rate*meanLoss*(tmax-Chara_material[id_mat].E0)/(Chara_material[id_mat].E0*tmax*log(w1));

	if(a1>nmaxCont)	
	{
		emean+=a1*e1;
		sig2e+=a1*e1*e1;
	}
	else
		if(a1>0.)
		{
			p1=(double)G4Poisson(a1);
			LossFluct+=p1*e1;
			if(p1>0.)
				LossFluct+=(1.-2.*Uniform())*e1;
		}
	if(a2>nmaxCont)
	{
		emean+=a2*e2;
		sig2e+=a2*e2*e2;
	}
	else
		if(a2>0.)
		{
			p2=(double)G4Poisson(a2);
			LossFluct+=p2*e2;
			if(p2>0.)
				LossFluct+=(1.-2.*Uniform())*e2;
		}
	if(a3>0.)
	{
		p3=a3;
		alfa=1.;
		if(a3>nmaxCont)
		{
			alfa=w1*(nmaxCont+a3)/(w1*nmaxCont+a3);
			alfa1=alfa*log(alfa)/(alfa-1.);
			nmean=a3*w1*(alfa-1.)/((w1-1.)*alfa);
			emean=emean+nmean*Chara_material[id_mat].E0*alfa1;
			sig2e=sig2e+Chara_material[id_mat].E0*Chara_material[id_mat].E0*nmean*(alfa-alfa1*alfa1);
			p3=a3-nmean;
		}
		w2=alfa*Chara_material[id_mat].E0;
		w=(tmax-w2)/tmax;
		nb=G4Poisson(p3);
		if(nb>0)
			for(i=0;i<nb;i++)	
				lossc+=w2/(1.-w*Uniform());
	}
	if(emean>0.)
	{
		sige=sqrt(sig2e);
		LossFluct+=Max(0.,Gaussian(emean,sige));
	}
	LossFluct+=lossc;

	return	LossFluct;
}

double eFluctuation(double meanLoss,double cutEnergy)
{
	int	nb,k;
	double	LossFluct=0.,lossc=0.,minLoss=10.*eV;
	double	rate=.55,fw=4.,nmaxCont=16.;
	double	tau,gamma,gamma2,beta2;
	double	F1,F2,E0,E1,E2,E1Log,E2Log,I,ILog;
	double	e1,e2,esmall,w,w1,w2,C,alfa,alfa1,namean;
	double	a1=0.,a2=0.,a3=0.,sa1;
	double	emean=0.,sig2e=0.,sige=0.,p1=0.,p2=0.,p3=0.;
	double	tmax=Min(cutEnergy,.5*Particle[id].Ekin);

	if(Particle[id].nature==positron)
		tmax=Min(cutEnergy,Particle[id].Ekin);

	if(meanLoss<minLoss)
		return	meanLoss;
	tau=Particle[id].Ekin/Particle[id].mass;
	gamma=tau+1.;
	gamma2=gamma*gamma;
	beta2=tau*(tau+2.)/gamma2;
	F1=Chara_material[id_mat].F1;
	F2=Chara_material[id_mat].F2;
	E0=Chara_material[id_mat].E0;
	E1=Chara_material[id_mat].E1;
	E2=Chara_material[id_mat].E2;
	E1Log=Chara_material[id_mat].E1Log;
	E2Log=Chara_material[id_mat].E2Log;
	I=Chara_material[id_mat].excitation_energy;
	ILog=Chara_material[id_mat].excitationLog;
	esmall=.5*sqrt(E0*I);  

	if(tmax<=E0)
		return	meanLoss;
	if(tmax>I) 
	{
		w2=log(2.*electron_mass_c2*beta2*gamma2)-beta2;
		if(w2>ILog)
		{
			C=meanLoss*(1.-rate)/(w2-ILog);
			a1=C*F1*(w2-E1Log)/E1;
			if(w2>E2Log)
				a2=C*F2*(w2-E2Log)/E2;
			if(a1<nmaxCont) 
			{ 
				sa1=sqrt(a1);
				if(Uniform()<exp(-sa1))
				{
					e1=esmall;
					a1=meanLoss*(1.-rate)/e1;
					a2=0.;
					e2=E2;
				} 
				else
				{
					a1=sa1 ;    
					e1=sa1*E1;
					e2=E2;
				}
			}
			else
			{
				a1/=fw;
				e1=fw*E1;
				e2=E2;
			}
		}   
	}
	w1=tmax/E0;
	if(tmax>E0) 
		a3=rate*meanLoss*(tmax-E0)/(E0*tmax*log(w1));
	if(a1>nmaxCont)
	{
		emean+=a1*e1;
		sig2e+=a1*e1*e1;
	}
	else
		if(a1>0.)
		{
			p1=(double)G4Poisson(a1);
			LossFluct+=p1*e1;
			if(p1>0.) 
				LossFluct+=(1.-2.*Uniform())*e1;
		}
	if(a2>nmaxCont)
	{
		emean+=a2*e2;
		sig2e+=a2*e2*e2;
	}
	else 
		if(a2>0.)
		{
			p2=(double)G4Poisson(a2);
			LossFluct+=p2*e2;
			if(p2>0.) 
				LossFluct+=(1.-2.*Uniform())*e2;
		}
	if(a3>0.)
	{
		p3=a3;
		alfa=1.;
		if(a3>nmaxCont)
		{
			alfa=w1*(nmaxCont+a3)/(w1*nmaxCont+a3);
			alfa1=alfa*log(alfa)/(alfa-1.);
			namean=a3*w1*(alfa-1.)/((w1-1.)*alfa);
			emean+=namean*E0*alfa1;
			sig2e+=E0*E0*namean*(alfa-alfa1*alfa1);
			p3=a3-namean;
		}	
		w2=alfa*E0;
		w=(tmax-w2)/tmax;
		nb=G4Poisson(p3);
		if(nb>0)
		for(k=0;k<nb;k++) 
			lossc+=w2/(1.-w*Uniform());
	}
	if(emean>0.)
	{
		sige=sqrt(sig2e);
		LossFluct+=Max(0.,Gaussian(emean,sige));
	}
	LossFluct+=lossc;

	return	LossFluct;
}

double hLoss(double LossLength,double cutEnergyElectron)
{
	double	perteTot=0.,perteIoni=0.;
	perteIoni=GetDedx(Particle[id].Ekin,1);
	perteTot=LossLength*perteIoni;
	if(perteTot>Particle[id].Ekin*xi)
		perteTot=LossApproximation(LossLength);
	if(bool.fluctuation==true)
		perteTot=hFluctuation(perteTot,cutEnergyElectron,LossLength);
	if(Particle[id].Ekin-perteTot<1.*eV)
	{
		perteTot=Particle[id].Ekin;
		bool.laststep=true;
	}	
	Particle[id].Ekin-=perteTot;
	return	perteTot;	
}

double eLoss(double LossLength,double cutEnergyElectron)
{
	double	perteTot=0.,perteBrem=0.,perteIoni=0.;
	perteIoni=GetDedx(Particle[id].Ekin,1);
	perteBrem=GetDedx(Particle[id].Ekin,2);
	perteTot=LossLength*(perteBrem+perteIoni);
	if(perteTot>Particle[id].Ekin*xi)
		perteTot=LossApproximation(LossLength);
	if(bool.fluctuation==true&&perteIoni>0.)	
		perteTot=eFluctuation(perteTot,cutEnergyElectron);
	if(Particle[id].Ekin-perteTot<1.*eV)
	{
		perteTot=Particle[id].Ekin;
		bool.laststep=true;
	}	
	Particle[id].Ekin-=perteTot;
	return	perteTot;	
}

double GlobalLoss(double LossLength,double cutEnergyElectron)
{
	int	particle;
	double	GlobalLoss=0.;

	particle=Particle[id].nature;
	if(bool.perte==false)
		return;
	switch(particle)
	{
		case electron:
		case positron:
			GlobalLoss=eLoss(LossLength,cutEnergyElectron);
		break;
		case proton:
			GlobalLoss=hLoss(LossLength,cutEnergyElectron);
		break;
	}
	if(bool.verbose==true)
		printf("Energie : %lf Step : %2.2E\n",Particle[id].Ekin,LossLength);

	return	GlobalLoss;
}

double hBetheIoniCrossSectionPerAtom(double CutEnergyElectron,int i)
{
	double	Sigma_I=0.,tmax,MaxEnergy;
	double	energytot2,beta2;
	
	tmax=Tmax();
	MaxEnergy=Min(tmax,MaxKinEnergy);
	if(MaxEnergy>CutEnergyElectron)
	{
		energytot2=Particle[id].Ekin+Particle[id].mass;
		energytot2=energytot2*energytot2;
		beta2=Particle[id].Ekin*(Particle[id].Ekin+2.*Particle[id].mass)/energytot2;
		Sigma_I=1./CutEnergyElectron-1./MaxEnergy-beta2*log(MaxEnergy/CutEnergyElectron)/tmax;
		if(Particle[id].Ekin>2.*MeV)
			Sigma_I+=.5*(MaxEnergy-CutEnergyElectron)/energytot2;
		Sigma_I*=twopi_mc2_rcl2*Particle[id].charge*Particle[id].charge/beta2;
	}
	Sigma_I*=Chara_material[id_mat].Comp_material[i].Z;
	return	Sigma_I;
}

double eIoniCrossSectionPerAtom(double cutEnergy,int i)
{
	double	Cross=0.;
	double	tmax=Min(1.*GeV,.5*Particle[id].Ekin);
	double	xmin,xmax,gamma,gamma2,beta2,g;
	double	y,y2,y12,b1,b2,y122,b4,b3;

	if(Particle[id].nature==positron)
		tmax=Min(1.*GeV,Particle[id].Ekin);

	if(cutEnergy<tmax) 
	{
		xmin=cutEnergy/Particle[id].Ekin;
		xmax=tmax/Particle[id].Ekin;
		gamma=Particle[id].Ekin/electron_mass_c2+1.;
		gamma2=gamma*gamma;
		beta2=1.-1./gamma2;
		if(Particle[id].nature==electron)
		{	
			g=(2.*gamma-1.)/gamma2;
			Cross=((xmax-xmin)*(1.-g+1./(xmin*xmax)+1./((1.-xmin)*(1.-xmax)))
					 -g*log(xmax*(1.-xmin)/(xmin*(1.-xmax))))/beta2;
		}
		else
		{
      y=1./(1.+gamma);
      y2=y*y;
      y12=1.-2.*y;
      b1=2.-y2;
      b2=y12*(3.+y2);
      y122=y12*y12;
      b4=y122*y12;
      b3=b4+y122;
      Cross=(xmax-xmin)*(1./(beta2*xmin*xmax)+b2-.5*b3*(xmin+xmax)
	    		 +b4*(xmin*xmin+xmin*xmax+xmax*xmax)/3.)-b1*log(xmax/xmin);
		}
		Cross*=twopi_mc2_rcl2/Particle[id].Ekin;
	}
	Cross*=Chara_material[id_mat].Comp_material[i].Z;
	return	Cross;
}

double hBetheIoniCrossSection(double CutEnergyElectron)
{
	int	i;
	double	Sigma_pv=0.;
	for(i=0;i<Chara_material[id_mat].Nbr_elmt;i++)
		Sigma_pv+=Chara_material[id_mat].Comp_material[i].D*Chara_material[id_mat].NbOfAtomsPerVolume
						*hBetheIoniCrossSectionPerAtom(CutEnergyElectron,i);
	return	Sigma_pv;
}

double eIoniCrossSection(double CutEnergyElectron)
{
	int	i;
	double	CrossTotale=0.;
	
	for(i=0;i<Chara_material[id_mat].Nbr_elmt;i++)
		CrossTotale+=Chara_material[id_mat].Comp_material[i].D*Chara_material[id_mat].NbOfAtomsPerVolume
								*eIoniCrossSectionPerAtom(CutEnergyElectron,i);
	return	CrossTotale;
}

double eBremCrossSectionPerAtom(double Z,double cut)
{
	int	i,j,iz=0,NZ=8,Nsig=11;
  double	Cross=0.;
  double	ksi=2.,alfa=1.;
  double	csigh=.127,csiglow=.25,asiglow=.02*MeV;
  double	Tlim=10.*MeV;
  double	xlim=1.2,delz=1.E6,absdelz;
  double	xx,fs;
  double	x,eta,alfap,K,a1,a3,a5;
	
	if (Particle[id].Ekin<1.*keV||Particle[id].Ekin<cut)
  	return	Cross;

	for(i=0;i<NZ;i++)
	{
		absdelz=fabs(Z-ZZ[i]); 
		if(absdelz<delz)
		{
			iz=i;
			delz=absdelz;
		}
	}

	xx=log10(Particle[id].Ekin);
	fs=1.;
	if(xx<=xlim)
	{
		fs=coefsig[iz][Nsig-1];
		for(j=Nsig-2;j>=0;j--)
			fs=fs*xx+coefsig[iz][j];
		if(fs<0.)
			fs=0.;
	}
	Cross=Z*(Z+ksi)*(1.-csigh*exp(log(Z)/4.))*pow(log(Particle[id].Ekin/cut),alfa);

	if(Particle[id].Ekin<=Tlim)
		Cross*=exp(csiglow*log(Tlim/Particle[id].Ekin))*(1.+asiglow/(sqrt(Z)*Particle[id].Ekin));

	if(Particle[id].nature==positron)
	{
		K=132.9416*eV;
		a1=4.15E-1;
		a3=2.1E-3;
		a5=54.E-5;
		x=log(Particle[id].Ekin/(K*Z*Z));
		eta=.5+atan(a1*x+a3*pow(x,3.)+a5*pow(x,5.))/pi;
		alfap=(1.-eta)/eta;
		Cross*=eta*pow((1.-cut/Particle[id].Ekin),alfap);;
	}

	Cross*=fs/N_avogadro;
	if(Cross<0.)
		Cross=0.;

	return	Cross;
}

double eBremCrossSectionPerVolume(double cutEnergy)
{
	int	i,n,nn,nmax=100;
	double	Cross=0.;
	double	kmax,kmin,vmin,vmax,totalEnergy,kp2;
	double	u,fac,c,v,dv,y;
	double	tmax=Min(MaxKinEnergy,Particle[id].Ekin);
	double	cut=Max(cutEnergy,.1*keV);
	double	fsig=0.;
	double	highKinEnergy=1.*GeV;
	double	probsup=1.;
	double	MigdalConstant=elec_radius*hbarc*hbarc*4.*pi/(electron_mass_c2*electron_mass_c2);
	
	if(cut>=tmax)
		return Cross;

	for(i=0;i<Chara_material[id_mat].Nbr_elmt;i++) 
	{
		Cross+=Chara_material[id_mat].Comp_material[i].D*Chara_material[id_mat].NbOfAtomsPerVolume
					*eBremCrossSectionPerAtom(Chara_material[id_mat].Comp_material[i].Z,cut);
		if(tmax<Particle[id].Ekin) 
			Cross-=Chara_material[id_mat].Comp_material[i].D*Chara_material[id_mat].NbOfAtomsPerVolume
						*eBremCrossSectionPerAtom(Chara_material[id_mat].Comp_material[i].Z,tmax);
	}

	kmax=tmax;
	kmin=cut;
	totalEnergy=Particle[id].Ekin+electron_mass_c2;
	kp2=MigdalConstant*totalEnergy*totalEnergy*Chara_material[id_mat].ElecDensity;
	vmin=log(kmin);
	vmax=log(kmax) ;
	nn=(int)(nmax*(vmax-vmin)/(log(highKinEnergy)-vmin));
	if(nn>0)
	{
		dv=(vmax-vmin)/nn;
		v=vmin-dv;
		for(n=0;n<=nn;n++)
		{
			v+=dv;  
			u=exp(v);              
			//fac=SupressionFunction(material,Particle[id].Ekin,u);		//LPM flag is off
			fac=1.;
			y=u/kmax;
			fac*=(4.-4.*y+3.*y*y)/3.;
			fac*=probsup*(u*u/(u*u+kp2))+1.-probsup;
			if((n==0)||(n==nn)) 
				c=.5;
			else
				c=1.;
			fac*=c;
			fsig+=fac;
		}
		y=kmin/kmax;
		fsig*=dv/(-4.*log(y)/3.-4.*(1.-y)/3.+0.5*(1.-y*y));
	} 
	else 
		fsig=1.;
	if(fsig>1.)
		fsig=1.;
	Cross*=fsig;

	return Cross;
}

double eBremCrossSection(double CutEnergyGamma)
{
	int	i;
	double CrossTotale=0.;
	CrossTotale=eBremCrossSectionPerVolume(CutEnergyGamma);
	return	CrossTotale;							
}

double pAnniCrossSection()
{
	double	Cross=0.;
	double	tau,gamma,gamma2,bg,bg2;
	
	tau=Particle[id].Ekin/electron_mass_c2;
	gamma=tau+1.;
	gamma2=gamma*gamma;
	bg2=tau*(tau+2.);
	bg=sqrt(bg2);
	Cross=((gamma2+4*gamma+1.)*log(gamma+bg)-(gamma+3.)*bg)/(bg2*(gamma+1.))
			 *Chara_material[id_mat].ElecDensity*pi*elec_radius*elec_radius;
	return Cross; 
}

double hMscCrossSectionPerAtom(int i)	
{
	int	iZ,iT;
	double	epsmin=1.E-4,epsmax=1.E10,Tlim=10.*MeV;	
	double	Tau,c,w,tau,electron_kinetic_energy,electron_total_energy;
	double	bg2,beta2,Z23,eps,bg2lim,beta2lim,Z1,Z2,ratZ,T,E,b2small,b2big,ratb2;
	double	c1,c2,cc1,cc2,correction,Sigma_p_a;
	                  	                 			 
	electron_kinetic_energy=Particle[id].Ekin;
	if(Particle[id].mass>electron_mass_c2)
	{
		Tau=Particle[id].Ekin/Particle[id].mass;
		c=Particle[id].Ekin*(Tau+2.)/(electron_mass_c2*(Tau+1.));
		w=c-2.;
		tau=.5*(w+sqrt(w*w+4.*c));
		electron_kinetic_energy=electron_mass_c2*tau;
	}			
	electron_total_energy=electron_kinetic_energy+electron_mass_c2;
	bg2=electron_kinetic_energy*(electron_total_energy+electron_mass_c2)
		/(electron_mass_c2*electron_mass_c2);
	beta2=electron_kinetic_energy*(electron_total_energy+electron_mass_c2)
			/(electron_total_energy*electron_total_energy);
  Z23=exp(2.*log(Chara_material[id_mat].Comp_material[i].Z)/3.);
  eps=2.*electron_mass_c2*electron_mass_c2*Bohr_radius*Bohr_radius/(hbarc*hbarc)*bg2/Z23;
  				
  if(eps<epsmin)		
  	Sigma_p_a=2.*eps*eps;
  else 
  	if(eps<epsmax)		
  		Sigma_p_a=log(1.+2.*eps)-2.*eps/(1.+2.*eps);
  	else
  		Sigma_p_a=log(2.*eps)-1.+1./eps;
  Sigma_p_a*=Chara_material[id_mat].Comp_material[i].Z*Chara_material[id_mat].Comp_material[i].Z/(beta2*bg2);

	iZ=14;
  while((iZ>=0)&&(Zdat[iZ]>=Chara_material[id_mat].Comp_material[i].Z))	
  	iZ-=1;
  if(iZ==14)
  	iZ=13;
  if(iZ==-1)
  	iZ=0;
	Z1=Zdat[iZ];
  Z2=Zdat[iZ+1];
  ratZ=(Chara_material[id_mat].Comp_material[i].Z-Z1)*(Chara_material[id_mat].Comp_material[i].Z+Z1)/((Z2-Z1)*(Z2+Z1));
  
  if(electron_kinetic_energy<=Tlim)
  {
  	iT=21;
  	while((iT>=0)&&(Tdat[iT]>=electron_kinetic_energy))	
  		iT-=1;
  	if(iT==21)																						
  		iT=20;
  	if(iT==-1)																						
  		iT=0;
  	T=Tdat[iT];
  	E=T+electron_mass_c2;
  	b2small=T*(E+electron_mass_c2)/(E*E);
  	T=Tdat[iT+1];
  	E=T+electron_mass_c2;
  	b2big=T*(E+electron_mass_c2)/(E*E);
	 	ratb2=(beta2-b2small)/(b2big-b2small);
  	
  	if(Particle[id].charge>0.)
  	{
  		c1=cpositron[iZ][iT];
  		c2=cpositron[iZ+1][iT];
  		cc1=c1+ratZ*(c2-c1);
  		c1=cpositron[iZ][iT+1];
  		c2=cpositron[iZ+1][iT+1];
  		cc2=c1+ratZ*(c2-c1);
  		correction=cc1+ratb2*(cc2-cc1);
  		Sigma_p_a*=2.*pi*elec_radius*elec_radius/correction;
  	}	
  }
  else 
  {
  	bg2lim=Tlim*(Tlim+2.*electron_mass_c2)/(electron_mass_c2*electron_mass_c2);
  	beta2lim=Tlim*(Tlim+2.*electron_mass_c2)/((Tlim+electron_mass_c2)*(Tlim+electron_mass_c2));
  	c1=bg2lim*sig0[iZ]*(1.+hecorr[iZ]*(beta2-beta2lim))/bg2;								
  	c2=bg2lim*sig0[iZ+1]*(1.+hecorr[iZ+1]*(beta2-beta2lim))/bg2;
  	if((Chara_material[id_mat].Comp_material[i].Z>=Z1)&&(Chara_material[id_mat].Comp_material[i].Z<=Z2))
  		Sigma_p_a=c1+ratZ*(c2-c1);
  	else
  		if(Chara_material[id_mat].Comp_material[i].Z<Z1)
  			Sigma_p_a=Chara_material[id_mat].Comp_material[i].Z*Chara_material[id_mat].Comp_material[i].Z*c1/(Z1*Z1);
  		else	
  			if(Chara_material[id_mat].Comp_material[i].Z>Z2)
  				Sigma_p_a=Chara_material[id_mat].Comp_material[i].Z*Chara_material[id_mat].Comp_material[i].Z*c2/(Z2*Z2);
  }
  return	Sigma_p_a;
}

double eMscCrossSectionPerAtom(double AtomicNumber)
{
	int	iZ=14,iT=21;
	double	Cross=0.;
	double	eKin,eTot,T,E;
	double	beta2,bg2,b2big,b2small,ratb2,Z23,tau,w;
	double	Z1,Z2,ratZ;
	double	c,c1,c2,cc1,cc2,corr;
	double	Tlim=10.*MeV;
	double	sigmafactor=2.*pi*elec_radius*elec_radius;
	double	epsfactor=2.*electron_mass_c2*electron_mass_c2*Bohr_radius*Bohr_radius/(hbarc*hbarc);
	double	eps,epsmin=1.e-4,epsmax=1.e10;
	double	beta2lim=Tlim*(Tlim+2.*electron_mass_c2)/((Tlim+electron_mass_c2)*(Tlim+electron_mass_c2));
	double	bg2lim=Tlim*(Tlim+2.*electron_mass_c2)/(electron_mass_c2*electron_mass_c2);
	double	ChargeSquare=Particle[id].charge*Particle[id].charge;
	
	Z23=2.*log(AtomicNumber)/3.; 
	Z23=exp(Z23);

	eKin=Particle[id].Ekin;
	if(Particle[id].mass>electron_mass_c2)
	{
		tau=Particle[id].Ekin/Particle[id].mass;
		c=Particle[id].mass*tau*(tau+2.)/(electron_mass_c2*(tau+1.));
		w=c-2.;
		tau=.5*(w+sqrt(w*w+4.*c));
		eKin=electron_mass_c2*tau;
	}

	eTot=eKin+electron_mass_c2;
	beta2=eKin*(eTot+electron_mass_c2)/(eTot*eTot);
	bg2=eKin*(eTot+electron_mass_c2)/(electron_mass_c2*electron_mass_c2);
	eps=epsfactor*bg2/Z23;
	if(eps<epsmin)
		Cross=2.*eps*eps;
	else 
		if(eps<epsmax)
			Cross=log(1.+2.*eps)-2.*eps/(1.+2.*eps);
		else
			Cross=log(2.*eps)-1.+1./eps;
	Cross*=ChargeSquare*AtomicNumber*AtomicNumber/(beta2*bg2);

	while((iZ>=0)&&(Zdat[iZ]>=AtomicNumber))
		iZ-=1;
	if(iZ==14)
		iZ=13;
	if(iZ==-1)
		iZ=0;
	Z1=Zdat[iZ];
	Z2=Zdat[iZ+1];
	ratZ=(AtomicNumber-Z1)*(AtomicNumber+Z1)/((Z2-Z1)*(Z2+Z1));

	if(eKin<=Tlim) 
	{
		while((iT>=0)&&(Tdat[iT]>=eKin))
			iT-=1;
		if(iT==21)
			iT=20;
		if(iT==-1)
			iT=0;
		T=Tdat[iT];
		E=T+electron_mass_c2;
		b2small=T*(E+electron_mass_c2)/(E*E);
		T=Tdat[iT+1]; 
		E=T+electron_mass_c2;
		b2big=T*(E+electron_mass_c2)/(E*E);
		ratb2=(beta2-b2small)/(b2big-b2small);

		if(Particle[id].charge<0.)
		{
			c1=celectron[iZ][iT];
			c2=celectron[iZ+1][iT];
			cc1=c1+ratZ*(c2-c1);
			c1=celectron[iZ][iT+1];
			c2=celectron[iZ+1][iT+1];
			cc2=c1+ratZ*(c2-c1);
			corr=cc1+ratb2*(cc2-cc1);
			Cross*=sigmafactor/corr;
		}
		else              
		{
			c1=cpositron[iZ][iT];
			c2=cpositron[iZ+1][iT];
			cc1=c1+ratZ*(c2-c1);
			c1=cpositron[iZ][iT+1];
			c2=cpositron[iZ+1][iT+1];
			cc2=c1+ratZ*(c2-c1);
			corr=cc1+ratb2*(cc2-cc1);
			Cross*=sigmafactor/corr;
		}
	}
	else
	{
		c1=bg2lim*sig0[iZ]*(1.+hecorr[iZ]*(beta2-beta2lim))/bg2;
		c2=bg2lim*sig0[iZ+1]*(1.+hecorr[iZ+1]*(beta2-beta2lim))/bg2;
		if((AtomicNumber>=Z1)&&(AtomicNumber<=Z2))
			Cross=c1+ratZ*(c2-c1);
		else 
			if(AtomicNumber<Z1)
				Cross=AtomicNumber*AtomicNumber*c1/(Z1*Z1);
			else 
				if(AtomicNumber>Z2)
					Cross=AtomicNumber*AtomicNumber*c2/(Z2*Z2);
	}
	return	Cross;
}

double hMscCrossSection()
{
	int i;
	double Sigma;
	Sigma=0.;
	for(i=0;i<Chara_material[id_mat].Nbr_elmt;i++)
	 Sigma+=Chara_material[id_mat].NbOfAtomsPerVolume*Chara_material[id_mat].Comp_material[i].D*hMscCrossSectionPerAtom(i);
	return	Sigma;
}

double eMscCrossSection()
{
	int	i;
	double CrossTotale=0.;
	
	for(i=0;i<Chara_material[id_mat].Nbr_elmt;i++)
		CrossTotale+=Chara_material[id_mat].Comp_material[i].D*Chara_material[id_mat].NbOfAtomsPerVolume
								*eMscCrossSectionPerAtom(Chara_material[id_mat].Comp_material[i].Z);
	
	return	CrossTotale;							
}

double gTransformToGeom(double TPath,double currentRange,double currentLambda,double currentEnergy,double *par1,double *par2)
{
	double	ZPath,zmean;
	double	tausmall=1.E-20,taulim=1.E-6,tlimitminfix=1.E-6*mm;
	double	dtrl=5./100.;
	double	tau,t1,lambda1;
	double	par3;
	
	*par1=-1.;
	*par2=par3=0.;
	ZPath=TPath;
	if(TPath<tlimitminfix)
		return	ZPath;
	if(TPath>currentRange)
		TPath=currentRange;
	tau=TPath/currentLambda;	
	if((tau<=tausmall)/*||insideskin*/) 
	{
		ZPath=TPath;
		if(ZPath>currentLambda) 
			ZPath=currentLambda;
		return	ZPath;
	}
  zmean=TPath;
	if(TPath<currentRange*dtrl) 
	{
		if(tau<taulim) 
			zmean=TPath*(1.-0.5*tau);
		else
			zmean=currentLambda*(1.-exp(-tau));
	} 
	else 
		if(currentEnergy<Particle[id].mass)  
		{
			*par1=1./currentRange;
			*par2=1./(*par1*currentLambda);
			par3=1.+*par2;
			if(TPath<currentRange)
				zmean=(1.-exp(par3*log(1.-TPath/currentRange)))/(*par1*par3);
			else
				zmean=1./(*par1*par3);
		} 
		else 
		{
			t1=GetEnergy(currentRange-TPath);
			lambda1=1./GetLambda(t1,1);
			*par1=(currentLambda-lambda1)/(currentLambda*TPath);
			*par2=1./(*par1*currentLambda);
			par3=1.+*par2;
			zmean=(1.-exp(par3*log(lambda1/currentLambda)))/(*par1*par3);
		}
	ZPath=zmean;
	if(ZPath>currentLambda) 
		ZPath=currentLambda;
	return	ZPath;
}

double eSimpleScattering(double xmeanth,double x2meanth)
{
  double	a=(2.*xmeanth+9.*x2meanth-3.)/(2.*xmeanth-3.*x2meanth+1.);
  double	prob=(a+2.)*xmeanth/a;
	double	cth=1.;
  if(Uniform()<prob)
    cth=-1.+2.*exp(log(Uniform())/(a+1.));
  else
    cth=-1.+2.*Uniform();
  return	cth;
}

double hCosineTheta(double trueStepLength,double currentRange,double currentLambda,double currentEnergy,double *currentTau,double par1,double par2)
{
	int	i,n_poisson;
	double	Cosinus_theta=1.;
	double	tau,mean,xmeanth,stepmin,taulim=1.E-6,taubig=8.,tausmall=1.E-20,dtrl=5./100.;
	double	theta0,c_highland,betacp,y,b,b1,bx,eb1,ebx,sth,c,c1,x0,f1x0,f2x0,
					a=1.,ea=0.,eaa=1.,xsi=3.,xmean1=1.,xmean2=0.,prob=0.,qprob=1.;
	double	tm,ascr,ascr1,bp1,bm1,ct,st,phi,sx=0.,sy=0.,sz=0.;	
	
	tau=trueStepLength/currentLambda;
	if(trueStepLength>=currentRange*dtrl)
	{
		if(par1*trueStepLength<1.)
			tau=-par2*log(1.-par1*trueStepLength);
		else
			if(1.-Particle[id].Ekin/currentEnergy>taulim)
				tau=taubig;	
	}
	*currentTau=tau;
	if(tau>=taubig)
		Cosinus_theta=-1.+2.*Uniform();
	else
		if(tau>=tausmall)
		{
			c_highland=13.6*MeV;
			betacp=sqrt(currentEnergy*(currentEnergy+2.*Particle[id].mass)*Particle[id].Ekin*(Particle[id].Ekin+2.*Particle[id].mass)
						/((currentEnergy+Particle[id].mass)*(Particle[id].Ekin+Particle[id].mass)));
			y=trueStepLength/Chara_material[id_mat].RadL;
			theta0=c_highland*sqrt(y)/betacp;
			y=log(y);				
			theta0*=sqrt(1.+y*(.105+.0035*y))*(1.-.24/(Chara_material[id_mat].Zeff*(Chara_material[id_mat].Zeff+1.)));
			if(theta0<tausmall)
				return	Cosinus_theta;
			sth=sin(.5*theta0);
			a=.25/(sth*sth);
			xmeanth=exp(-tau);
			c=3.;
			c1=c-1.;
			x0=1.-xsi/a;
			if(x0<0.)
			{
				b=exp(tau);
				bx=b-1.;
				b1=b+1.;
				ebx=exp((c1)*log(bx));
				eb1=exp((c1)*log(b1));
			}
			else
			{
				c=2.4-.027*exp(2.*log(Chara_material[id_mat].Zeff)/3.);
				if(c==2.)
					c=2.+taulim;
				if(c<=1.)
					c=1.+taulim;
				c1=c-1.;
				ea=exp(-xsi);
				eaa=1.-ea;
				xmean1=1.-(1.-(1.+xsi)*ea)/(eaa*a);
				b=1.+(c-xsi)/a;
				b1=b+1.;
				bx=c/a;
				ebx=exp((c1)*log(bx));
				eb1=exp((c1)*log(b1));
				xmean2=(x0*eb1+ebx-(eb1*bx-b1*ebx)/(c-2.))/(eb1-ebx);
				f1x0=a*ea/eaa;
				f2x0=c1*eb1*ebx/(eb1-ebx)/exp(c*log(bx));
				prob=f2x0/(f1x0+f2x0);
				qprob=(f1x0+f2x0)*xmeanth/(f2x0*xmean1+f1x0*xmean2);
			}
			
			if(Uniform()<qprob)
			{
				if(Uniform()<prob)
					Cosinus_theta=1.+log(ea+Uniform()*eaa)/a;
				else
					Cosinus_theta=b-b1*bx/exp(log(ebx-Uniform()*(ebx-eb1))/c1);
			}
			else
				Cosinus_theta=-1.+2.*Uniform();							
		}

	return Cosinus_theta;		
}

double eCosineTheta(double trueStep,double currentRange,double currentLambda,double currentEnergy,double *currentTau,double par1,double par2)
{
	double	costh,sinth;
	double	tau,taubig=8.,tausmall=1.E-16,taulim=1.E-6;
	double	c,c1,x0,b,bx,b1,ebx,eb1;
	double	prob=0.,qprob=1.;
	double	a=1.,ea=0.,eaa=1.;
	double	xmeanth,xmean1=1.,xmean2=0.,x2meanth;
	double	dtrl=5./100.;
	double	xsi=3.;
  double	theta0,theta0max=pi/6.,y,corr,betacp,c_highland=13.6*MeV;
  double	f1x0,f2x0;

	costh=1.;
	tau=trueStep/currentLambda;
	if(trueStep>=currentRange*dtrl) 
	{
		if((par1*trueStep)<1.)
			tau=-par2*log(1.-par1*trueStep);
		else 
			if((1.-Particle[id].Ekin/currentEnergy)>taulim)
				tau=taubig;
	}
	*currentTau=tau;
	if(tau>=taubig) 
		costh=-1.+2.*Uniform();
	else 
		if(tau>=tausmall)
		{
			x0=1.;	b=2.;	b1=3.;
			bx=1.;eb1=3.;ebx=1.;
			prob=1.;	qprob=1.;
			xmeanth=exp(-tau);
			x2meanth=(1.+2.*exp(-2.5*tau))/3.;
			if(1.-Particle[id].Ekin/currentEnergy>.5)
			{
				costh=eSimpleScattering(xmeanth,x2meanth);
				return	costh;
			}	

		betacp=sqrt(currentEnergy*(currentEnergy+2.*Particle[id].mass)
					*Particle[id].Ekin*(Particle[id].Ekin+2.*Particle[id].mass)
					/((currentEnergy+Particle[id].mass)*(Particle[id].Ekin+Particle[id].mass)));
		y=trueStep/Chara_material[id_mat].RadL;
		theta0=c_highland*fabs(Particle[id].charge)*sqrt(y)/betacp;
		y=log(y);
		corr=(1.-8.778E-2/Chara_material[id_mat].Zeff)*(.87+.03*log(Chara_material[id_mat].Zeff))
				+(4.078E-2+1.7315E-4*Chara_material[id_mat].Zeff)*(.87+.03*log(Chara_material[id_mat].Zeff))*y;                
		theta0*=corr ;                                               
		if(theta0*theta0<tausmall) 
			return	costh;
		if(theta0>theta0max)
		{
			costh=eSimpleScattering(xmeanth,x2meanth);
			return	costh;
		}	
	
		sinth=sin(.5*theta0);
		a=.25/(sinth*sinth);
		ea=exp(-xsi);
		eaa=1.-ea ;
		xmean1=1.-(1.-(1.+xsi)*ea)/(a*eaa);
		x0=1.-xsi/a;
		if(xmean1<=.999*xmeanth)
		{
			costh=eSimpleScattering(xmeanth,x2meanth);
			return	costh;
		}	

		c=2.943-.197*log(Chara_material[id_mat].Zeff+1.)
		+(.0987-.0143*log(Chara_material[id_mat].Zeff+1.))*y;

		if(fabs(c-3.)<.001)  
			c=3.001;      
		if(fabs(c-2.)<.001)  
			c=2.001;      
		if(fabs(c-1.)<.001)  
			c=1.001;      
		c1=c-1.;

		b=1.+(c-xsi)/a;
		b1=b+1.;
		bx=c/a;
		eb1=exp(c1*log(b1));
		ebx=exp(c1*log(bx));
		xmean2=(x0*eb1+ebx-(eb1*bx-b1*ebx)/(c-2.))/(eb1-ebx);
		f1x0=a*ea/eaa;
		f2x0=c1*eb1/(bx*(eb1-ebx));
		prob=f2x0/(f1x0+f2x0);
		qprob=xmeanth/(prob*xmean1+(1.-prob)*xmean2);
		if(Uniform()<qprob)
		{
			if(Uniform()<prob)
				costh=1.+log(ea+Uniform()*eaa)/a;
			else
				costh=b-b1*bx/exp(log(ebx+(eb1-ebx)*Uniform())/c1);
		}
		else
			costh=-1.+2.*Uniform();
	} 
	return	costh;
}

double gGeomLengthLimit(double gPath,double cStep,double currentLambda,double currentRange,double par1,double par3)
{
	double	tPath;
	double	tausmall=1.E-16;
	
	par3=1.+par3;
	tPath=gPath;
	if(gPath>currentLambda*tausmall)
	{
		if(par1<0.)
			tPath=-currentLambda*log(1.-gPath/currentLambda);
		else 
		{
			if(par1*par3*gPath<1.)
				tPath=(1.-exp(log(1.-par1*par3*gPath)/par3))/par1;
			else 
				tPath=currentRange;
		}
	}
	if(tPath<gPath)
		tPath=gPath;
	return	tPath;
}

void gLatCorrection(struct Vector currentDir,double tPath,double zPath,double currentTau,double phi,double sinth)
{
	double	safety,latcorr,etau,rmean,rmax,Phi,psi,lambdaeff;
	double	kappa=2.5,tlimitminfix=1.E-6*mm,taulim=1.E-6,tausmall=1.E-16,taubig=8.,geomMin=1.E-6*mm;
	struct	Vector	latDir;
	
	lambdaeff=tPath/currentTau;
	safety=GetSafety(Particle[id].Position);
	if(bool.verbose==true)
		printf("Safety : %lf\n",safety);
	if(bool.lateral==true&&safety>tlimitminfix) 
	{
		rmean=0.;
		if((currentTau>=tausmall)/* && !insideskin*/) 
		{
			if(currentTau<taulim) 
				rmean=kappa*pow(currentTau,3.)*(1.-(kappa+1.)*currentTau*.25)/6.;
			else 
			{
				etau=0.;
				if(currentTau<taubig) 
					etau=exp(-currentTau);
				rmean=-kappa*currentTau;
				rmean=-exp(rmean)/(kappa*(kappa-1.));
				rmean+=currentTau-(kappa+1.)/kappa+kappa*etau/(kappa-1.);
			}
			if(rmean>0.) 
				rmean=2.*lambdaeff*sqrt(rmean/3.);
			else
				rmean=0.;
		}
		rmax=(tPath-zPath)*(tPath+zPath);
		if(rmax<0.)
			rmax=0.;
		else
			rmax=sqrt(rmax);
		if(rmean>=rmax) 
			rmean=rmax;

		if(rmean<=geomMin)
			return;
		if(rmean>0.)
		{
			if((currentTau>=tausmall) /*&& !insideskin*/)
			{
				if(currentTau<taulim)
					latcorr=lambdaeff*kappa*currentTau*currentTau*(1.-(kappa+1.)*currentTau/3.)/3.;
				else
				{
					etau=0.;
					if(currentTau<taubig) 
						etau=exp(-currentTau);
					latcorr=-kappa*currentTau;
					latcorr=exp(latcorr)/(kappa-1.);
					latcorr+=1.-kappa*etau/(kappa-1.);
					latcorr*=2.*lambdaeff/3. ;
				}
			}	
			if(latcorr>rmean) 
				latcorr=rmean;
			Phi=0.;
			if(fabs(rmean*sinth)<latcorr)
				Phi=2.*pi*Uniform();
			else
			{
				psi=acos(latcorr/(rmean*sinth));
				if(Uniform()<.5)
					Phi=phi+psi;
				else
					Phi=phi-psi;
			}
			latDir.X=cos(Phi);
			latDir.Y=sin(Phi);
			latDir.Z=0.;
			latDir=RotateUz(currentDir,latDir);
			if(rmean>safety)
				rmean=safety*.99;
			Particle[id].Position.X+=latDir.X*rmean;
			Particle[id].Position.Y+=latDir.Y*rmean;
			Particle[id].Position.Z+=latDir.Z*rmean;
		}
	}
}

void hMscScattering(double tPath,double zPath,double currentRange,double currentLambda,double currentEnergy,double par1,double par2)
{
	double	costh,sinth,phi,dirx,diry,currentTau;
	struct Vector Dir,currentDir;
	
	costh=hCosineTheta(tPath,currentRange,currentLambda,currentEnergy,&currentTau,par1,par2);
	if(fabs(costh>1.))
		return;	
	sinth=sqrt((1.-costh)*(1.+costh));
	phi=2.*pi*Uniform();
	dirx=sinth*cos(phi);
	diry=sinth*sin(phi);
	Dir.X=dirx;	
	Dir.Y=diry;	
	Dir.Z=costh;
	Particle[id].Position.X+=Particle[id].Momentum.X*zPath;
	Particle[id].Position.Y+=Particle[id].Momentum.Y*zPath;
	Particle[id].Position.Z+=Particle[id].Momentum.Z*zPath;
	currentDir=Particle[id].Momentum;
	Particle[id].Momentum=RotateUz(Particle[id].Momentum,Dir);
	if(bool.lateral==true)
		gLatCorrection(currentDir,tPath,zPath,currentTau,phi,sinth);
}

void eMscScattering(double tPath,double zPath,double currentRange,double currentLambda,double currentEnergy,double par1,double par2)
{
	double	costh,sinth,phi,currentTau;
	double	tlimitminfix=1.E-10*mm,taulim=1.E-6,tausmall=1.E-16;
	struct	Vector	Dir,currentDir;

	if((Particle[id].Ekin<0.)||(tPath<=tlimitminfix)||(tPath/tausmall<currentLambda)) 
		return;

	costh=eCosineTheta(tPath,currentRange,currentLambda,currentEnergy,&currentTau,par1,par2);
	if(fabs(costh)>1.) 
		return;
	if(costh<1.-1000.*tPath/currentLambda&&Particle[id].Ekin>20.*MeV)
	{ 
		do{
			costh=1.+2.*log(Uniform())*tPath/currentLambda;
		}while(costh<-1.);
	}
	sinth=sqrt((1.-costh)*(1.+costh));
	phi=2.*pi*Uniform();
	Dir.X=sinth*cos(phi);
	Dir.Y=sinth*sin(phi);
	Dir.Z=costh;
	Particle[id].Position.X+=Particle[id].Momentum.X*zPath;
	Particle[id].Position.Y+=Particle[id].Momentum.Y*zPath;
	Particle[id].Position.Z+=Particle[id].Momentum.Z*zPath;
	currentDir=Particle[id].Momentum;
	Particle[id].Momentum=RotateUz(Particle[id].Momentum,Dir);
	if(bool.lateral==true)
		gLatCorrection(currentDir,tPath,zPath,currentTau,phi,sinth);
}

double GlobalMscScattering(double GeomPath,double CutStep,double CurrentRange,double CurrentLambda,double CurrentEnergy,double CutEnergyElectron,double par1,double par2)
{
	int	particle;
	double	Dose,TruePath,zPath,tausmall=1.E-16;
	particle=Particle[id].nature;
	if(bool.multiplescattering==false)
	{
		if(GeomPath<CutStep)
		{
			Dose=GlobalLoss(GeomPath,CutEnergyElectron);
			MapDose(false,Dose,Particle[id].Position);
			RendDose(false,Dose,Particle[id].Position);
		}
		Particle[id].Position.X+=Particle[id].Momentum.X*GeomPath;
		Particle[id].Position.Y+=Particle[id].Momentum.Y*GeomPath;
		Particle[id].Position.Z+=Particle[id].Momentum.Z*GeomPath;
		return	GeomPath;
	}
	
	if(GeomPath==CutStep)
		zPath=gTransformToGeom(GeomPath,CurrentRange,CurrentLambda,CurrentEnergy,&par1,&par2);
	else
	{
		zPath=GeomPath;
		TruePath=gGeomLengthLimit(GeomPath,CutStep,CurrentLambda,CurrentRange,par1,par2);
		GeomPath=TruePath;
		Dose=GlobalLoss(TruePath,CutEnergyElectron);
		MapDose(false,Dose,Particle[id].Position);
		RendDose(false,Dose,Particle[id].Position);
	}	
	if(bool.laststep==false)
		switch(particle)
		{
			case electron:
			case positron:
				eMscScattering(GeomPath,zPath,CurrentRange,CurrentLambda,CurrentEnergy,par1,par2);
			break;
			case proton:
				hMscScattering(GeomPath,zPath,CurrentRange,CurrentLambda,CurrentEnergy,par1,par2);
			break;
		}
	else
	{
		Particle[id].Position.X+=Particle[id].Momentum.X*zPath;
		Particle[id].Position.Y+=Particle[id].Momentum.Y*zPath;
		Particle[id].Position.Z+=Particle[id].Momentum.Z*zPath;
	}
	if(bool.verbose==true)
	{
		printf("zPath : %2.2E, TruePath : %2.2E, GeomPath : %2.2E\n",zPath,TruePath,GeomPath);
		printf("Position : ");
		PrintVector(Particle[id].Position);
		printf("Moment : ");
		PrintVector(Particle[id].Momentum);
	}
	return	TruePath;
}

void hSampleSecondarieElectron(double CutEnergy)
{
	double	deltaEnergy;
	double	tmax,tmin,q,g,x,x1,x2,eps,beta2,f,f1,fmax,totMom,deltaMom,cth,sth,phi;
	struct	Vector	ElecDir;

	eps=.8426*GeV;
	tmax=Min(Tmax(),MaxKinEnergy);
	tmin=Max(Chara_material[id_mat].excitation_energy,CutEnergy);
	beta2=Particle[id].Ekin*(Particle[id].Ekin+2.*Particle[id].mass)/pow((Particle[id].Ekin+Particle[id].mass),2.);
	f=0.;
	f1=0.;
	fmax=1.+.5*tmax*tmax/pow((Particle[id].Ekin+Particle[id].mass),2.);
	do
	{
		q=Uniform();
		deltaEnergy=tmin*tmax/(tmin*(1.-q)+tmax*q);
		if(Particle[id].Ekin>2.*MeV)
			f1=.5*deltaEnergy*deltaEnergy/pow((Particle[id].Ekin+Particle[id].mass),2.);
		f=f1+1.-beta2*deltaEnergy/tmax;
	}while(fmax*Uniform()>f);
	x=2.*electron_mass_c2/(eps*eps)*deltaEnergy;
	if(x>1.E-6&&Particle[id].Ekin>2.*MeV)
	{
		x1=x+1.;
		g=1./(x1*x1);
		x2=.5*electron_mass_c2*deltaEnergy/(Particle[id].mass*Particle[id].mass);
		g=g*(1.+(Magmom_proton*Magmom_proton-1.)*(x2-f1/f)/(1.+x2));
		if(g<Uniform())	
			return;
	}
	totMom=(Particle[id].Ekin+Particle[id].mass)*sqrt(beta2);
	deltaMom=sqrt(deltaEnergy*(deltaEnergy+2.*electron_mass_c2));
	cth=deltaEnergy*(Particle[id].Ekin+Particle[id].mass+electron_mass_c2)/(deltaMom*totMom);
	if(fabs(cth>1.))
		cth=1.;
	sth=sqrt((1.-cth)*(1.+cth));
	phi=2.*pi*Uniform();

	ElecDir.X=sth*cos(phi);
	ElecDir.Y=sth*sin(phi);
	ElecDir.Z=cth;
	ElecDir=RotateUz(Particle[id].Momentum,ElecDir);
	Particle[id].Momentum=CorrUnit(Particle[id].Momentum,ElecDir,totMom,deltaMom);
	if(deltaEnergy>Particle[id].Ekin)
		deltaEnergy=Particle[id].Ekin;
	Particle[id].Ekin-=deltaEnergy;	

	if(bool.secondaries==true&&(bool.tertiaries==true||id==0))
	{
		if(idsec<MAX_PART)
		{
			GenerateNewParticle(electron,deltaEnergy,ElecDir);
			// fprintf(DeltaRay,"%8.8lf %E %E %E %E %E %E\n",
			// Particle[idsec].Ekin,Particle[idsec].Position.X,Particle[idsec].Position.Y,Particle[idsec].Position.Z,
			// Particle[idsec].Momentum.X,Particle[idsec].Momentum.Y,Particle[idsec].Momentum.Z);	
		}
		else
		{
			printf("## Overflow of secondaries ##\n");
			bool.secondaries=false;
			bool.tertiaries=false;
		}	
	}
	else
		MapDose(false,deltaEnergy,Particle[id].Position);			
}

void eSampleSecondarieElectron(double CutEnergy)
{
	double	totalEnergy,deltaEnergy,totMom,deltaMom;
	double	xmin,xmax,gamma,gamma2,beta2;
	double	x,z,q,grej,g,y;
	double	y2,y12,y122,b1,b2,b3,b4;
	double	cost,sint,phi;
	double	tmax=Min(1.*GeV,.5*Particle[id].Ekin);
	double	tmin=CutEnergy;
	struct	Vector	ElecDir;

	if(Particle[id].nature==positron)
		tmax=Min(1.*GeV,Particle[id].Ekin);
	if(tmin>=tmax) 
		return;

	totalEnergy=Particle[id].Ekin+electron_mass_c2;
	totMom=sqrt(Particle[id].Ekin*(totalEnergy+ electron_mass_c2));
	xmin=tmin/Particle[id].Ekin;
	xmax=tmax/Particle[id].Ekin;
	gamma=totalEnergy/electron_mass_c2;
	gamma2=gamma*gamma;
	beta2=1.-1./gamma2;
	if(Particle[id].nature==electron)
	{	
		g=(2.*gamma-1.)/gamma2;
		y=1.-xmax;
		grej=1.-g*xmax+xmax*xmax*(1.-g+(1.-g*y)/(y*y));
		do{
			q=Uniform();
			x=xmin*xmax/(xmin*(1.-q)+xmax*q);
			y=1.-x;
			z=1.-g*x+x*x*(1.-g+(1.-g*y)/(y*y));
		}while(grej*Uniform()>z);
	}
	else
	{
    y=1./(1.+gamma);
    y2=y*y;
    y12=1.-2.*y;
    b1=2.-y2;
    b2=y12*(3.+y2);
    y122=y12*y12;
    b4=y122*y12;
    b3=b4+y122;
    y=xmax*xmax;
    grej=1.+(y*y*b4-xmin*xmin*xmin*b3+y*b2-xmin*b1)*beta2; 
    do{
    	q=Uniform();
    	x=xmin*xmax/(xmin*(1.-q)+xmax*q);
    	y=x*x;
    	z=1.+(y*y*b4-x*y*b3+y*b2-x*b1)*beta2; 
    }while(grej*Uniform()>z);
	}

	deltaEnergy=x*Particle[id].Ekin;
	deltaMom=sqrt(deltaEnergy*(deltaEnergy+2.*electron_mass_c2));
	cost=deltaEnergy*(totalEnergy+electron_mass_c2)/(deltaMom*totMom);
	sint=1.-cost*cost;
	if(sint>0.) 
		sint=sqrt(sint);
	phi=2.*pi*Uniform();
	
	ElecDir.X=sint*cos(phi);
	ElecDir.Y=sint*sin(phi);
	ElecDir.Z=cost;
	ElecDir=RotateUz(Particle[id].Momentum,ElecDir);
	Particle[id].Ekin-=deltaEnergy;
	if(Particle[id].Ekin>DBL_MIN) 
		Particle[id].Momentum=CorrUnit(Particle[id].Momentum,ElecDir,totMom,deltaMom);
	if(bool.secondaries==true&&(bool.tertiaries==true||id==0))
	{
		if(idsec<MAX_PART)
		{
			GenerateNewParticle(electron,deltaEnergy,ElecDir);
			// fprintf(DeltaRay,"%8.8lf %E %E %E %E %E %E\n",
			// Particle[idsec].Ekin,Particle[idsec].Position.X,Particle[idsec].Position.Y,Particle[idsec].Position.Z,
			// Particle[idsec].Momentum.X,Particle[idsec].Momentum.Y,Particle[idsec].Momentum.Z);	
		}
		else
		{
			printf("## Overflow of secondaries ##\n");
			bool.secondaries=false;
			bool.tertiaries=false;
		}
	}
	else
	 	MapDose(false,deltaEnergy,Particle[id].Position);		
}

void GlobalSampleSecondarieElectron(double cutEnergyElectron)
{
	int	particle;
	particle=Particle[id].nature;
	switch(particle)
	{
		case electron:
		case positron:
			eSampleSecondarieElectron(cutEnergyElectron);
		break;
		case proton:
			hSampleSecondarieElectron(cutEnergyElectron);
		break;
	}	
}

double ScreenFunction1(double ScreenVariable)
{
	double screenVal;
	if(ScreenVariable>1.)
		screenVal=42.24-8.368*log(ScreenVariable+.952);
	else
		screenVal=42.392-ScreenVariable*(7.796-1.961*ScreenVariable);
	return	screenVal;
} 

double ScreenFunction2(double ScreenVariable)
{
	double screenVal;
	if(ScreenVariable>1.)
		screenVal=42.24-8.368*log(ScreenVariable+.952);
	else
		screenVal=41.734-ScreenVariable*(6.484-1.25*ScreenVariable);
	return	screenVal;
} 

double RejectionFunction(double value,double rej1,double rej2,double rej3,double ratio,double z)
{
	double	argument=(1.+value)*(1.+value);
	return	(4.+log(rej3+(z/argument)))*((4.*ratio*value/argument)-rej1)+rej2;
}

double AngleDistribution(double initial_energy,double final_energy,double Z)
{
	double	Theta=0.;
	double	initialTotalEnergy=(initial_energy+electron_mass_c2)/electron_mass_c2;
	double	finalTotalEnergy=(final_energy+electron_mass_c2)/electron_mass_c2;
	double	EnergyRatio=finalTotalEnergy/initialTotalEnergy;
	double	gMaxEnergy=(pi*initialTotalEnergy)*(pi*initialTotalEnergy);
	double	z,rejection_argument1,rejection_argument2,rejection_argument3;
	double	gfunction0,gfunction1,gfunctionEmax,gMaximum;	
	double	rand,gfunctionTest,randTest;

	z=.00008116224*(pow(Z,1./3.)+pow(Z+1,1./3.));
	rejection_argument1=(1.+EnergyRatio*EnergyRatio); 
	rejection_argument2=-2.*EnergyRatio+3.*rejection_argument1;
	rejection_argument3=((1.-EnergyRatio)/(2.*initialTotalEnergy*EnergyRatio))*((1.-EnergyRatio)/(2.*initialTotalEnergy*EnergyRatio));
	gfunction0=RejectionFunction(0.,rejection_argument1,rejection_argument2,rejection_argument3,EnergyRatio,z);
	gfunction1=RejectionFunction(1.,rejection_argument1,rejection_argument2,rejection_argument3,EnergyRatio,z);
	gfunctionEmax=RejectionFunction(gMaxEnergy,rejection_argument1,rejection_argument2,rejection_argument3,EnergyRatio,z);
	gMaximum=Max(gfunction0,gfunction1);
	gMaximum=Max(gMaximum,gfunctionEmax);

	do{
		rand=Uniform();
		rand/=(1.-rand+1./gMaxEnergy);
		gfunctionTest=RejectionFunction(rand,rejection_argument1,rejection_argument2,rejection_argument3,EnergyRatio,z);
		randTest=Uniform();
	}while(randTest*gMaximum>gfunctionTest);
	Theta=sqrt(rand)/initialTotalEnergy;
	
	return	Theta;
}

int RandomAtom(double CutEnergyGamma)
{
	int i,indice,tmp;
	double U,rval;
	tmp=Chara_material[id_mat].Nbr_elmt-1;
	rval=Uniform()*Chara_material[id_mat].Comp_material[tmp].D*eBremCrossSectionPerAtom(Chara_material[id_mat].Comp_material[tmp].Z,CutEnergyGamma);
	for(i=0;i<Chara_material[id_mat].Nbr_elmt;i++)
	{
		U=Chara_material[id_mat].Comp_material[i].D*eBremCrossSectionPerAtom(Chara_material[id_mat].Comp_material[i].Z,CutEnergyGamma);
		if(rval<=U)
		{
			indice=i;
			break;
		}
	}
	return	indice;
}	
	
void eSampleSecondarieGamma(double cutEnergy)
{
	int	ind;
	double	gammaEnergy,totalEnergy;
	double	xmin,xmax,kappa,epsilmin,epsilmax;
	double	lnZ,FZ,Z3,ZZ,F1,F2,theta,sint,phi;
	double	tmin=cutEnergy;
	double	tmax=Min(MaxKinEnergy,Particle[id].Ekin);
	double	MigdalFactor,MigdalConstant=elec_radius*hbarc*hbarc*4.*pi/(electron_mass_c2*electron_mass_c2);
	double	x,xm,epsil,greject,migdal,grejmax,q,U,U2;
	double	ah,bh,screenvar,screenmin,screenfac=0.;
	double
		ah10= 4.67733E+00,ah11=-6.19012E-01,ah12= 2.02225E-02,
		ah20=-7.34101E+00,ah21= 1.00462E+00,ah22=-3.20985E-02,
		ah30= 2.93119E+00,ah31=-4.03761E-01,ah32= 1.25153E-02;
	double
		bh10= 4.23071E+00,bh11=-6.10995E-01,bh12= 1.95531E-02,
		bh20=-7.12527E+00,bh21= 9.69160E-01,bh22=-2.74255E-02,
		bh30= 2.69925E+00,bh31=-3.63283E-01,bh32= 9.55316E-03;
	double
		al00=-2.05398E+00,al01= 2.38815E-02,al02= 5.25483E-04,
		al10=-7.69748E-02,al11=-6.91499E-02,al12= 2.22453E-03,
		al20= 4.06463E-02,al21=-1.01281E-02,al22= 3.40919E-04;
	double
		bl00= 1.04133E+00,bl01=-9.43291E-03,bl02=-4.54758E-04,
		bl10= 1.19253E-01,bl11= 4.07467E-02,bl12=-1.30718E-03,
		bl20=-1.59391E-02,bl21= 7.27752E-03,bl22=-1.94405E-04;
	double	ah1,ah2,ah3,bh1,bh2,bh3;
	double	al1,al2,al0,bl1,bl2,bl0;
	double	tlow = 1.*MeV;
	double	totMom;
	struct	Vector	GamDir;

	if(tmin>=tmax)
		return;
	ind=RandomAtom(cutEnergy);
	
	Z3=pow(Chara_material[id_mat].Comp_material[ind].Z,1./3.);
	lnZ=3.*log(Z3);
	FZ=lnZ*(4.-.55*lnZ);
	ZZ=pow(Chara_material[id_mat].Comp_material[ind].Z*(Chara_material[id_mat].Comp_material[ind].Z+1.),1./3.);

	totalEnergy=Particle[id].Ekin+electron_mass_c2;
	xmin=tmin/Particle[id].Ekin;
	xmax=tmax/Particle[id].Ekin;
	kappa=0.;
	if(xmax>=1.)
		xmax=1.;
	else
		kappa=log(xmax)/log(xmin);
	epsilmin=tmin/totalEnergy;
	epsilmax=tmax/totalEnergy;
	MigdalFactor=Chara_material[id_mat].ElecDensity*MigdalConstant/(epsilmax*epsilmax);
	U=log(Particle[id].Ekin/electron_mass_c2);
	U2=U*U;

	if(Particle[id].Ekin>tlow)
	{
		ah1=ah10+ZZ*(ah11+ZZ*ah12);
		ah2=ah20+ZZ*(ah21+ZZ*ah22);
		ah3=ah30+ZZ*(ah31+ZZ*ah32);
		bh1=bh10+ZZ*(bh11+ZZ*bh12);
		bh2=bh20+ZZ*(bh21+ZZ*bh22);
		bh3=bh30+ZZ*(bh31+ZZ*bh32);
		ah=1.+(ah1*U2+ah2*U+ah3)/(U2*U);
		bh=.75+(bh1*U2+bh2*U+bh3)/(U2*U);
		screenfac=136.*electron_mass_c2/(Z3*totalEnergy);
		screenmin=screenfac*epsilmin/(1.-epsilmin);
		F1=Max(ScreenFunction1(screenmin)-FZ,0.);
		F2=Max(ScreenFunction2(screenmin)-FZ,0.);
		grejmax=(F1-epsilmin*(F1*ah-bh*epsilmin*F2))/(42.392-FZ);
	} 
	else 
	{  
		al0=al00+ZZ*(al01+ZZ*al02);
		al1=al10+ZZ*(al11+ZZ*al12);
		al2=al20+ZZ*(al21+ZZ*al22);
		bl0=bl00+ZZ*(bl01+ZZ*bl02);
		bl1=bl10+ZZ*(bl11+ZZ*bl12);
		bl2=bl20+ZZ*(bl21+ZZ*bl22);
		ah=al0+al1*U+al2*U2;
		bh=bl0+bl1*U+bl2*U2;
		grejmax=Max(1.+xmin*(ah+bh*xmin),1.+ah+bh);
		xm=-ah/(2.*bh);
		if(xmin<xm&&xm<xmax) 
			grejmax=Max(grejmax,1.+xm*(ah+bh*xm));
	}

	if(Particle[id].Ekin>tlow) 
	{
		do{
			q=Uniform();
			x=pow(xmin,q+kappa*(1.-q));
			epsil=x*Particle[id].Ekin/totalEnergy;
			screenvar=screenfac*epsil/(1.-epsil);
			F1=Max(ScreenFunction1(screenvar)-FZ,0.);
			F2=Max(ScreenFunction2(screenvar)-FZ,0.);
			migdal=(1.+MigdalFactor)/(1.+MigdalFactor/(x*x));
			greject=migdal*(F1-epsil*(ah*F1-bh*epsil*F2))/(42.392-FZ);
		}while(greject<Uniform()*grejmax);
	} 
	else 
	{  
		do{
			q=Uniform();
			x=pow(xmin,q+kappa*(1.-q));
			migdal=(1.+MigdalFactor)/(1.+MigdalFactor/(x*x));  
			greject=migdal*(1.+x*(ah+bh*x));
		}while(greject<Uniform()*grejmax);
	}
	gammaEnergy=x*Particle[id].Ekin; 

	theta=AngleDistribution(totalEnergy,totalEnergy-gammaEnergy,Chara_material[id_mat].Comp_material[ind].Z);
	sint=sin(theta);
	phi=2.*pi*Uniform();
	GamDir.X=sint*cos(phi);
	GamDir.Y=sint*sin(phi);
	GamDir.Z=cos(theta);
	GamDir=RotateUz(Particle[id].Momentum,GamDir);
	totMom=sqrt(Particle[id].Ekin*(totalEnergy+electron_mass_c2));
	Particle[id].Momentum=CorrUnit(Particle[id].Momentum,GamDir,totMom,gammaEnergy);
	Particle[id].Ekin=Particle[id].Ekin-gammaEnergy;
	MapDose(false,gammaEnergy,Particle[id].Position);		
}

int CorrectionMat(int indice)
{
	if(indice==22||indice==75||indice==98)
		indice=G4_Adipose_Tissue;
	else
	if(indice==0||indice==15||indice==16||indice==29||indice==104)
		indice=G4_Air;
	else
	if(indice==100)
		indice=G4_B100_Bone;
	else
	if(indice==23||indice==84)
		indice=G4_Blood;
	else
	if(indice==5||indice==70)
		indice=G4_Bone_Compact;
	else
	if(indice==77||indice==85||indice==91||indice==96||indice==101||indice==105||indice==108||indice==111||indice==114||indice==118||
		 indice==83||indice==89||indice==95||indice==97||indice==103||indice==107||indice==109||indice==112||indice==117||indice==120||indice==124)
		indice=G4_Brain;
	else
	if(indice==121)
		indice=G4_Eye_Lens;
	else
	if(indice==200)
		indice=G4_MS20_Tissue;
	else
	if(indice==9||indice==78)
		indice=G4_Muscle_Skeletal;
	else
	if(indice==110||indice==116||indice==119)
		indice=G4_Muscle_Sucrose;
	else
	if(indice==1||indice==30||indice==82||indice==113)
		indice=G4_Skin;
	else
	if(indice==2||indice==3||indice==4||indice==26||indice==72||indice==74||indice==92||indice==106||indice==115||indice==122||indice==123||indice==125)
		indice=G4_Water;
	return	indice;
}

void PhantomReader()
{
	char	char_mat;
	int	i,j,k,tmp1,tmp2,tmp3,tmp_mat;
	int	nbVoxel,idVoxel;
	double	dtmp1,dtmp2,dtmp3;
	
	for(i=0;i<MAX_MAT;i++)
		Mater[i]=false;


	// Phantom_file=fopen("Phantom/atlas_380x200x208/atlas_208x200x380.dim","r");
	// fscanf(Phantom_file,"%d %d %d\n",&tmp1,&tmp2,&tmp3);
	// Voxel.bin.x=tmp1;	Voxel.bin.y=tmp2;	Voxel.bin.z=tmp3;
	// Voxel.Dimension.X=.1*mm;	Voxel.Dimension.Y=.1*mm;	Voxel.Dimension.Z=.1*mm;
	// nbVoxel=Voxel.bin.x*Voxel.bin.y*Voxel.bin.z;
	// Vecteur_Voxel_mat=(int *)malloc(sizeof(int)*nbVoxel);
	// for(i=0;i<nbVoxel;i++)
	// 	Vecteur_Voxel_mat[i]=0;
	// fclose(Phantom_file);

	// Phantom_file=fopen("Phantom/atlas_380x200x208/atlas_208x200x380.ima","rb");
	// for(k=0;k<Voxel.bin.z;k++)	for(j=0;j<Voxel.bin.y;j++)	for(i=0;i<Voxel.bin.x;i++)
	// {	
	// 	fread(&char_mat,sizeof(char),1,Phantom_file);
	// 	tmp_mat=(int)char_mat;
	// 	if(tmp_mat==4||tmp_mat==5||tmp_mat==6||tmp_mat==7||tmp_mat==8||tmp_mat==12||tmp_mat==13||tmp_mat==15||tmp_mat==20)
	// 		tmp_mat=G4_Water;
	// 	idVoxel=k*Voxel.bin.x*Voxel.bin.y+j*Voxel.bin.x+i;
	// 	Vecteur_Voxel_mat[idVoxel]=tmp_mat;
	// 	Mater[tmp_mat]=true;
	// }
	// fclose(Phantom_file);

	//Phantom_file=fopen("Phantom/phantom.in","r");
	Phantom_file=fopen("Phantom/Phantom_water_cube.dat","r");
	//Phantom_file=fopen("Phantom/phantom_rubiks.in","r");
	//Phantom_file=fopen("Phantom/phantom_mush.in","r");
	//Phantom_file=fopen("Phantom/phantom_bone.in","r");
	//Phantom_file=fopen("Phantom/ascii_phantom_carcinoma_base.dat","r");
	//Phantom_file=fopen("Phantom/ascii_cube_.dat","r");
	//Phantom_file=fopen("Phantom/ascii_cube_eau.dat","r");
	if(!Phantom_file)
	{
		printf("/!\\ Error, no phantom file\n");
		exit(0);
	}
	fscanf(Phantom_file,"%d %d %d\n",&tmp1,&tmp2,&tmp3);
	Voxel.bin.x=tmp1;	Voxel.bin.y=tmp2;	Voxel.bin.z=tmp3;
	fscanf(Phantom_file,"%lf %lf %lf\n",&dtmp1,&dtmp2,&dtmp3);
	Voxel.Dimension.X=dtmp1*mm;	Voxel.Dimension.Y=dtmp2*mm;	Voxel.Dimension.Z=dtmp3*mm;
	nbVoxel=Voxel.bin.x*Voxel.bin.y*Voxel.bin.z;
	Vecteur_Voxel_mat=(int *)malloc(sizeof(int)*nbVoxel);
	for(i=0;i<nbVoxel;i++)
		Vecteur_Voxel_mat[i]=0;

	// for(k=0;k<20;k++)	
	// 	for(j=0;j<Voxel.bin.y;j++)	
	// 		for(i=0;i<Voxel.bin.x;i++)
	// 		{	
	// 			tmp_mat=1;
	// 			idVoxel=k*Voxel.bin.x*Voxel.bin.y+j*Voxel.bin.x+i;
	// 			Vecteur_Voxel_mat[idVoxel]=tmp_mat;
	// 			Mater[tmp_mat]=true;
	// 		}
	// for(k=20;k<25;k++)	
	// 	for(j=0;j<Voxel.bin.y;j++)	
	// 		for(i=0;i<Voxel.bin.x;i++)
	// 		{	
	// 			tmp_mat=2;
	// 			idVoxel=k*Voxel.bin.x*Voxel.bin.y+j*Voxel.bin.x+i;
	// 			Vecteur_Voxel_mat[idVoxel]=tmp_mat;
	// 			Mater[tmp_mat]=true;
	// 		}
	// for(k=25;k<Voxel.bin.z;k++)	
	// 	for(j=0;j<Voxel.bin.y;j++)	
	// 		for(i=0;i<Voxel.bin.x;i++)
	// 		{	
	// 			tmp_mat=6;
	// 			idVoxel=k*Voxel.bin.x*Voxel.bin.y+j*Voxel.bin.x+i;
	// 			Vecteur_Voxel_mat[idVoxel]=tmp_mat;
	// 			Mater[tmp_mat]=true;
	// 		}
	
	for(k=0;k<Voxel.bin.z;k++)	
		for(j=0;j<Voxel.bin.y;j++)	
			for(i=0;i<Voxel.bin.x;i++)
			{	
				fscanf(Phantom_file,"%d",&tmp_mat);
				tmp_mat=CorrectionMat(tmp_mat);
				idVoxel=k*Voxel.bin.x*Voxel.bin.y+j*Voxel.bin.x+i;
				Vecteur_Voxel_mat[idVoxel]=tmp_mat;
				Mater[tmp_mat]=true;
			}

	// fscanf(Phantom_file,"%d",&tmp_mat);
	// for(k=0;k<Voxel.bin.z;k++)	for(j=0;j<Voxel.bin.y;j++)	for(i=0;i<Voxel.bin.x;i++)
	// {	
	// 	idVoxel=k*Voxel.bin.x*Voxel.bin.y+j*Voxel.bin.x+i;
	// 	Vecteur_Voxel_mat[idVoxel]=tmp_mat;
	// 	Mater[tmp_mat]=true;
	// }
	fclose(Phantom_file);
}

void Energy_table()
{
	int	i;
	double	constant,slope,x,energy;
	
	constant=MinKinEnergy;
	slope=log(MaxKinEnergy/MinKinEnergy);
	for(i=0;i<Bins;i++)
	{
		x=(double)i;
		x/=(Bins-1);
		Energy_tab[id_mat][i]=constant*exp(slope*x)*MeV;
		hInvRange_tab[id_mat][i]=Energy_tab[id_mat][i];
		eInvRange_tab[id_mat][i]=Energy_tab[id_mat][i];
		pInvRange_tab[id_mat][i]=Energy_tab[id_mat][i];
	}
	max_ind=Bins;
}

void Range_table()
{
	int	i,j,n;
	double	energy,de,hsum,esum,psum,hDXDE=0.,eDXDE=0.,pDXDE=0.;
	
	i=0;
	n=100;
	hDXDE=hIoni_DEDX_tab[id_mat][i];
	if(hDXDE>0.)
		hDXDE=2.*Energy_tab[id_mat][i]/hDXDE;
	hRange_tab[id_mat][i]=hDXDE;		
	eDXDE=eIoni_DEDX_tab[id_mat][i]+eBrem_DEDX_tab[id_mat][i];
	if(eDXDE>0.)
		eDXDE=2.*Energy_tab[id_mat][i]/eDXDE;
	eRange_tab[id_mat][i]=eDXDE;		
	pDXDE=pIoni_DEDX_tab[id_mat][i]+pBrem_DEDX_tab[id_mat][i];
	if(pDXDE>0.)
		pDXDE=2.*Energy_tab[id_mat][i]/pDXDE;
	pRange_tab[id_mat][i]=pDXDE;		
	for(i=1;i<max_ind;i++)
	{
		de=(Energy_tab[id_mat][i]-Energy_tab[id_mat][i-1])/n;
		energy=Energy_tab[id_mat][i]+de*.5;
		hsum=0.;
		esum=0.;
		psum=0.;
		for(j=0;j<n;j++)
		{
			energy-=de;
			Particle[id].nature=proton;
			hDXDE=GetDedx(energy,1);
			if(hDXDE>0.)
				hsum+=de/hDXDE;
			Particle[id].nature=electron;
			eDXDE=GetDedx(energy,1)+GetDedx(energy,2);
			if(eDXDE>0.)
				esum+=de/eDXDE;
			Particle[id].nature=positron;
			pDXDE=GetDedx(energy,1)+GetDedx(energy,2);
			if(pDXDE>0.)
				psum+=de/pDXDE;
		}		
		hRange_tab[id_mat][i]=hRange_tab[id_mat][i-1]+hsum;
		eRange_tab[id_mat][i]=eRange_tab[id_mat][i-1]+esum;
		pRange_tab[id_mat][i]=pRange_tab[id_mat][i-1]+psum;
	}
}

void hMSC_CrossSection_table()
{
	int	i;
	Particle[id].nature=proton;
	Particle[id].mass=m_proton;
	Particle[id].charge=1.;
	for(i=0;i<max_ind;i++)
	{
		Particle[id].Ekin=Energy_tab[id_mat][i];
		hMSC_CS_tab[id_mat][i]=hMscCrossSection();
	}
}

void eMSC_CrossSection_table()
{
	int	i;
	Particle[id].nature=electron;
	Particle[id].mass=electron_mass_c2;
	Particle[id].charge=-1.;
	for(i=0;i<max_ind;i++)
	{
		Particle[id].Ekin=Energy_tab[id_mat][i];
		eMSC_CS_tab[id_mat][i]=eMscCrossSection();
	}
}

void pMSC_CrossSection_table()
{
	int	i;
	Particle[id].nature=positron;
	Particle[id].mass=electron_mass_c2;
	Particle[id].charge=1.;
	for(i=0;i<max_ind;i++)
	{
		Particle[id].Ekin=Energy_tab[id_mat][i];
		pMSC_CS_tab[id_mat][i]=eMscCrossSection();
	}
}

void hIoni_CrossSection_table(double cutEnergyElectron)
{
	int	i;
	Particle[id].nature=proton;
	Particle[id].mass=m_proton;
	Particle[id].charge=1.;
	for(i=0;i<max_ind;i++)
	{
		Particle[id].Ekin=Energy_tab[id_mat][i];
		hIoni_CS_tab[id_mat][i]=hBetheIoniCrossSection(cutEnergyElectron);
	}
}

void eIoni_CrossSection_table(double cutEnergyElectron)
{
	int	i;
	Particle[id].nature=electron;
	Particle[id].mass=electron_mass_c2;
	Particle[id].charge=-1.;
	for(i=0;i<max_ind;i++)
	{
		Particle[id].Ekin=Energy_tab[id_mat][i];
		eIoni_CS_tab[id_mat][i]=eIoniCrossSection(cutEnergyElectron);
	}
}

void pIoni_CrossSection_table(double cutEnergyElectron)
{
	int	i;
	Particle[id].nature=positron;
	Particle[id].mass=electron_mass_c2;
	Particle[id].charge=1.;
	for(i=0;i<max_ind;i++)
	{
		Particle[id].Ekin=Energy_tab[id_mat][i];
		pIoni_CS_tab[id_mat][i]=eIoniCrossSection(cutEnergyElectron);
	}
}

void eBrem_CrossSection_table(double cutEnergyGamma)
{
	int	i;
	Particle[id].nature=electron;
	Particle[id].mass=electron_mass_c2;
	Particle[id].charge=-1.;
	for(i=0;i<max_ind;i++)
	{
		Particle[id].Ekin=Energy_tab[id_mat][i];
		eBrem_CS_tab[id_mat][i]=eBremCrossSection(cutEnergyGamma)*mm2;	//G4 internal unit;
	}
}

void pBrem_CrossSection_table(double cutEnergyGamma)
{
	int	i;
	Particle[id].nature=positron;
	Particle[id].mass=electron_mass_c2;
	Particle[id].charge=1.;
	for(i=0;i<max_ind;i++)
	{
		Particle[id].Ekin=Energy_tab[id_mat][i];
		pBrem_CS_tab[id_mat][i]=eBremCrossSection(cutEnergyGamma)*mm2;	//G4 internal unit;
	}
}

void pAnni_CrossSection_table()
{
	int	i;
	Particle[id].nature=positron;
	Particle[id].mass=electron_mass_c2;
	Particle[id].charge=1.;
	for(i=0;i<max_ind;i++)
	{
		Particle[id].Ekin=Energy_tab[id_mat][i];
		pAnni_CS_tab[id_mat][i]=pAnniCrossSection();
	}
}

void hIoni_DEDX_table(double cutEnergyElectron)
{
	int	i;
	double	Tlim=2.*MeV,SlTlim,ShTlim,ParamB;
	Particle[id].nature=proton;
	Particle[id].mass=m_proton;
	Particle[id].charge=1.;
	Particle[id].Ekin=Tlim;
	SlTlim=hBraggDEDX(cutEnergyElectron);
	ShTlim=hBetheIoniDEDX(cutEnergyElectron);
	ParamB=SlTlim/ShTlim-1.;
	for(i=0;i<max_ind;i++)
	{
		Particle[id].Ekin=Energy_tab[id_mat][i];
		if(Particle[id].Ekin<=Tlim)
			hIoni_DEDX_tab[id_mat][i]=hBraggDEDX(cutEnergyElectron);
		else
			hIoni_DEDX_tab[id_mat][i]=hBetheIoniDEDX(cutEnergyElectron)*(1.+ParamB*Tlim/Energy_tab[id_mat][i]);
	}
}

void eIoni_DEDX_table(double cutEnergyElectron)
{
	int	i;
	Particle[id].nature=electron;
	Particle[id].mass=electron_mass_c2;
	Particle[id].charge=-1.;
	for(i=0;i<max_ind;i++)
	{
		Particle[id].Ekin=Energy_tab[id_mat][i];
		eIoni_DEDX_tab[id_mat][i]=eIoniDEDX(cutEnergyElectron);
	}
}

void pIoni_DEDX_table(double cutEnergyElectron)
{
	int	i;
	Particle[id].nature=positron;
	Particle[id].mass=electron_mass_c2;
	Particle[id].charge=1.;
	for(i=0;i<max_ind;i++)
	{
		Particle[id].Ekin=Energy_tab[id_mat][i];
		pIoni_DEDX_tab[id_mat][i]=eIoniDEDX(cutEnergyElectron);
	}
}

void eBrem_DEDX_table(double cutEnergyGamma)
{
	int	i;
	Particle[id].nature=electron;
	Particle[id].mass=electron_mass_c2;
	Particle[id].charge=-1.;
	for(i=0;i<max_ind;i++)
	{
		Particle[id].Ekin=Energy_tab[id_mat][i];
		eBrem_DEDX_tab[id_mat][i]=eBremDEDX(cutEnergyGamma)*mm2;	//G4 internal unit
	}
}

void pBrem_DEDX_table(double cutEnergyGamma)
{
	int	i;
	Particle[id].nature=positron;
	Particle[id].mass=electron_mass_c2;
	Particle[id].charge=1.;
	for(i=0;i<max_ind;i++)
	{
		Particle[id].Ekin=Energy_tab[id_mat][i];
		pBrem_DEDX_tab[id_mat][i]=eBremDEDX(cutEnergyGamma)*mm2;	//G4 internal unit
	}
}

void PrintVector(struct Vector Vecteur)
{
	printf("(%8.8E;%8.8E;%8.8E)\n",Vecteur.X,Vecteur.Y,Vecteur.Z);
}

void PrintTable(double cutEnergyGamma,double cutEnergyElectron)
{
	int	i;
	if(bool.eperte==true||bool.ebrem==true||bool.emultiplescattering==true)
	{
		printf("# electron in material: %d\nMSC: %d, cut energy electron: %5.5E MeV, cut energy gamma: %5.5E MeV #\n",
						id_mat,bool.emultiplescattering,cutEnergyElectron,cutEnergyGamma);
		printf("	Energy			Dedx			 Range			CS Ioni			CS Brem			CS MSC\n");
		for(i=0;i<max_ind;i++)
			printf("%3.3E	%6.6E	%6.6E	%6.6E	%6.6E	%6.6E\n",
							Energy_tab[id_mat][i],eIoni_DEDX_tab[id_mat][i]+eBrem_DEDX_tab[id_mat][i],eRange_tab[id_mat][i],
							eIoni_CS_tab[id_mat][i],eBrem_CS_tab[id_mat][i],eMSC_CS_tab[id_mat][i]);
		printf("\n");
 	}
	if(bool.eperte==true||bool.ebrem==true||bool.emultiplescattering==true||bool.panni==true)
	{
		printf("# positron in material: %d\nMSC: %d, cut energy electron: %5.5E MeV, cut energy gamma: %5.5E MeV #\n",
						id_mat,bool.emultiplescattering,cutEnergyElectron,cutEnergyGamma);
		printf("	Energy			Dedx			 Range			CS Ioni			CS Brem			CS MSC			CS Anni\n");
		for(i=0;i<max_ind;i++)
			printf("%3.3E	%6.6E	%6.6E	%6.6E	%6.6E	%6.6E	%6.6E\n",
							Energy_tab[id_mat][i],pIoni_DEDX_tab[id_mat][i]+pBrem_DEDX_tab[id_mat][i],pRange_tab[id_mat][i],
							pIoni_CS_tab[id_mat][i],pBrem_CS_tab[id_mat][i],pMSC_CS_tab[id_mat][i],pAnni_CS_tab[id_mat][i]);
		printf("\n");
 	}
 	if(bool.hperte==true||bool.hmultiplescattering==true)
 	{
 		printf("# proton in material: %d\nMSC: %d, cut energy electron: %5.5E MeV #\n",
 						id_mat,bool.hmultiplescattering,cutEnergyElectron);
 		printf("	Energy			Dedx			 Range			CS Ioni			CS MSC\n");
 		for(i=0;i<max_ind;i++)
 			printf("%3.3E	%6.6E	%6.6E	%6.6E	%6.6E\n",
 							Energy_tab[id_mat][i],hIoni_DEDX_tab[id_mat][i],hRange_tab[id_mat][i],
 							hIoni_CS_tab[id_mat][i],hMSC_CS_tab[id_mat][i]);
 		printf("\n");
	}
}

void PrintGeometry()
{
	double	theta=Source.Angle;
	double	coord_x1,coord_y1,coord_z1,coord_x2,coord_y2,coord_z2;
	FILE *Target_geometry,*Detector_geometry,*Source_geometry;
	
	Target_geometry=fopen("../../Python/Volumes/Target.fit","w");
	coord_x1=-Target.Dimension.X;
	coord_y1=-Target.Dimension.Y;
	coord_z1=-Target.Dimension.Z;
	coord_x2=-coord_x1;
	coord_y2=-coord_y1;
	coord_z2=-coord_z1;
	fprintf(Target_geometry,"%lf %lf %lf\n",coord_x1,coord_y1,coord_z1); //0
	fprintf(Target_geometry,"%lf %lf %lf\n",coord_x1,coord_y1,coord_z2); //1
	fprintf(Target_geometry,"%lf %lf %lf\n",coord_x1,coord_y2,coord_z2); //2
	fprintf(Target_geometry,"%lf %lf %lf\n",coord_x1,coord_y2,coord_z1); //3
	fprintf(Target_geometry,"%lf %lf %lf\n",coord_x1,coord_y1,coord_z1); //0
	fprintf(Target_geometry,"%lf %lf %lf\n",coord_x2,coord_y1,coord_z1); //4
	fprintf(Target_geometry,"%lf %lf %lf\n",coord_x2,coord_y1,coord_z2); //5
	fprintf(Target_geometry,"%lf %lf %lf\n",coord_x1,coord_y1,coord_z2); //1
	fprintf(Target_geometry,"%lf %lf %lf\n",coord_x1,coord_y2,coord_z2); //2
	fprintf(Target_geometry,"%lf %lf %lf\n",coord_x2,coord_y2,coord_z2); //6
	fprintf(Target_geometry,"%lf %lf %lf\n",coord_x2,coord_y2,coord_z1); //7
	fprintf(Target_geometry,"%lf %lf %lf\n",coord_x2,coord_y1,coord_z1); //4
	fprintf(Target_geometry,"%lf %lf %lf\n",coord_x2,coord_y1,coord_z2); //5
	fprintf(Target_geometry,"%lf %lf %lf\n",coord_x2,coord_y2,coord_z2); //6
	fprintf(Target_geometry,"%lf %lf %lf\n",coord_x2,coord_y2,coord_z1); //7
	fprintf(Target_geometry,"%lf %lf %lf\n",coord_x1,coord_y2,coord_z1); //3
	fclose(Target_geometry);
	
	Detector_geometry=fopen("../../Python/Volumes/Detector.fit","w");
	coord_x1=Detecteur.Dimension.Z*sin(theta)-Detecteur.Dimension.X*cos(theta);
	coord_y1=-Detecteur.Dimension.Y;
	coord_z1=Detecteur.Dimension.Z*cos(theta)+Detecteur.Dimension.X*sin(theta);
	coord_x2=Detecteur.Dimension.Z*sin(theta)+Detecteur.Dimension.X*cos(theta);
	coord_y2=-coord_y1;
	coord_z2=Detecteur.Dimension.Z*cos(theta)-Detecteur.Dimension.X*sin(theta);
	fprintf(Detector_geometry,"%lf %lf %lf\n",coord_x1,coord_y1,coord_z1); //0
	fprintf(Detector_geometry,"%lf %lf %lf\n",coord_x1,coord_y2,coord_z1); //1
	fprintf(Detector_geometry,"%lf %lf %lf\n",coord_x2,coord_y2,coord_z2); //2
	fprintf(Detector_geometry,"%lf %lf %lf\n",coord_x2,coord_y1,coord_z2); //3
	fprintf(Detector_geometry,"%lf %lf %lf\n",coord_x1,coord_y1,coord_z1); //0
	fclose(Detector_geometry);
	
	Source_geometry=fopen("../../Python/Volumes/Source.fit","w");
	coord_x1=Source.Dimension.Z*sin(theta)-Source.Dimension.X*cos(theta);
	coord_y1=-Source.Dimension.Y;
	coord_z1=Source.Dimension.Z*cos(theta)+Source.Dimension.X*sin(theta);
	coord_x2=Source.Dimension.Z*sin(theta)+Source.Dimension.X*cos(theta);
	coord_y2=-coord_y1;
	coord_z2=Source.Dimension.Z*cos(theta)-Source.Dimension.X*sin(theta);
	fprintf(Source_geometry,"%lf %lf %lf\n",coord_x1,coord_y1,coord_z1); //0
	fprintf(Source_geometry,"%lf %lf %lf\n",coord_x1,coord_y2,coord_z1); //1
	fprintf(Source_geometry,"%lf %lf %lf\n",coord_x2,coord_y2,coord_z2); //2
	fprintf(Source_geometry,"%lf %lf %lf\n",coord_x2,coord_y1,coord_z2); //3
	fprintf(Source_geometry,"%lf %lf %lf\n",coord_x1,coord_y1,coord_z1); //0
	fclose(Source_geometry);
}

void Generate_table(double cutEnergyGamma,double cutEnergyElectron)
{
	int	i;
	for(i=0;i<MAX_RNG;i++)
	{
		Energy_tab[id_mat][i]=0.;
		hRange_tab[id_mat][i]=0.;
		hInvRange_tab[id_mat][i]=0.;
		eRange_tab[id_mat][i]=0.;
		eInvRange_tab[id_mat][i]=0.;
		pRange_tab[id_mat][i]=0.;
		pInvRange_tab[id_mat][i]=0.;
		hMSC_CS_tab[id_mat][i]=0.;
		eMSC_CS_tab[id_mat][i]=0.;
		pMSC_CS_tab[id_mat][i]=0.;
		hIoni_CS_tab[id_mat][i]=0.;
		eIoni_CS_tab[id_mat][i]=0.;
		pIoni_CS_tab[id_mat][i]=0.;
		eBrem_CS_tab[id_mat][i]=0.;
		pBrem_CS_tab[id_mat][i]=0.;
		hIoni_DEDX_tab[id_mat][i]=0.;
		eIoni_DEDX_tab[id_mat][i]=0.;
		pAnni_CS_tab[id_mat][i]=0.;
	}
	id=0;
	Energy_table();
	if(bool.hperte==true)
	{
		hIoni_DEDX_table(cutEnergyElectron);
		hIoni_CrossSection_table(cutEnergyElectron);
	}
	if(bool.eperte==true)
	{
		eIoni_DEDX_table(cutEnergyElectron);
		eIoni_CrossSection_table(cutEnergyElectron);
		pIoni_DEDX_table(cutEnergyElectron);
		pIoni_CrossSection_table(cutEnergyElectron);
	}	
	if(bool.ebrem==true)
	{
		eBrem_DEDX_table(cutEnergyGamma);
		eBrem_CrossSection_table(cutEnergyGamma);
		pBrem_DEDX_table(cutEnergyGamma);
		pBrem_CrossSection_table(cutEnergyGamma);
	}	
	if(bool.hperte==true||bool.eperte==true||bool.ebrem==true)
		Range_table();
	if(bool.hmultiplescattering==true)
		hMSC_CrossSection_table();							 
	if(bool.emultiplescattering==true)
	{	
		eMSC_CrossSection_table();
		pMSC_CrossSection_table();
	}
	if(bool.panni==true)
		pAnni_CrossSection_table();
}

int DefineProcess(int flag)
{
	int i,prc;

	Processus[0].CrossSection=GetLambda(Particle[id].Ekin,2);
	Processus[1].CrossSection=GetLambda(Particle[id].Ekin,3);
	if(Particle[id].nature==positron)
		Processus[2].CrossSection=GetLambda(Particle[id].Ekin,4);
	else
		Processus[2].CrossSection=0.;
	prc=0;
	
	switch(flag)
	{
		case true:
			Processus[0].PreviousCS=Processus[0].CrossSection;
			for(i=0;i<MAX_PRC;i++)
			{
				if(Processus[i].CrossSection==0.)
					Processus[i].TrueLength=DBL_MAX;
				else
					Processus[i].TrueLength=-log(Uniform())/Processus[i].CrossSection;
				if(Processus[prc].TrueLength>Processus[i].TrueLength)
					prc=i;													 
			}
		break;
		
		case false:
			for(i=0;i<MAX_PRC;i++)
			{
				if(Processus[i].CrossSection==0.)
					Processus[i].TrueLength=DBL_MAX;
				else
					if(Processus[i].PreviousCS==0.)
						Processus[i].TrueLength=-log(Uniform())/Processus[i].CrossSection;
					else	
						Processus[i].TrueLength*=Processus[i].PreviousCS/Processus[i].CrossSection;
				if(Processus[prc].TrueLength>Processus[i].TrueLength)
					prc=i;													 
				Processus[i].PreviousCS=Processus[i].CrossSection;
			}
		break;
	}
	return	prc;	
}

int Navigate(double *totalLength,double *freeLength,double cutEnergyElectron,double cutEnergyGamma,double EkinLimit,int *face_voxel)
{
	int	bool_loop,bool_cut,face_cube,prc,significant_loss;
	double	alpha,rho,cutStep,par1,par2;
	double	Dose;
	double	alongStepLength,lengthtoVertex,trueStepLength,trueGeomLength,safety,fragment;
	double	currentEnergy,currentRange,currentLambda;
	
	alpha=0.;
	rho=0.*mm;
	bool_loop=true;
	bool_cut=true;
	face_cube=0;
	cutStep=1.*mm;
	alongStepLength=*freeLength;
	bool.laststep=false;
	if(alongStepLength>0.)
		bool_cut=false;

	if(bool.verbose==true)
		printf("Materiau : %d\n",id_mat);

	if(bool.multiplescattering==false&&bool.perte==false)
	{
		safety=GetStepSide(Particle[id].Position,cutStep,&face_cube,false);
		RayTracing(safety);
		*freeLength=alongStepLength+safety;
		*totalLength+=safety;
		*face_voxel=face_cube;
		if(bool.ray_tracing==true&&Particle[id].level==0)
				fprintf(Ray_tracing_file,"%lf %lf %lf %d\n",Particle[id].Position.X,Particle[id].Position.Y,Particle[id].Position.Z,Particle[id].nature);
		return	false;
	}
	
	do{
		prc=DefineProcess(bool_cut);
		lengthtoVertex=VertexLength(Processus[prc].TrueLength,alongStepLength);
		currentRange=GetRange(Particle[id].Ekin);	currentEnergy=Particle[id].Ekin;	currentLambda=1./GetLambda(currentEnergy,1);
		cutStep=StepFunction(currentRange,alpha,rho);
		safety=GetStepSide(Particle[id].Position,cutStep,&face_cube,false);
		if(bool.verbose==true)
			printf("currentRange : %2.2E cutStep : %2.2E lengthtoVertex : %2.2E safety : %2.2E\n",currentRange,cutStep,lengthtoVertex,safety);
		if(lengthtoVertex>cutStep)
		{
			significant_loss=true;
			trueStepLength=cutStep;
		}
		else
		{
			significant_loss=false;
			trueStepLength=lengthtoVertex;
		}
		if(trueStepLength>safety)
		{
			if(bool.multiplescattering==false)
				bool_loop=false;
			else
			{
				trueGeomLength=gTransformToGeom(trueStepLength,currentRange,currentLambda,currentEnergy,&par1,&par2);
				if(trueGeomLength>safety)
					bool_loop=false;
			}
		}
		else
			if(bool.verbose==true)
				printf("Step : %2.2E\n",trueStepLength);
		if(bool_loop==true)
		{
			switch(significant_loss)
			{
				case true:
					Dose=GlobalLoss(trueStepLength,cutEnergyElectron);
					GlobalMscScattering(trueStepLength,cutStep,currentRange,currentLambda,currentEnergy,cutEnergyElectron,par1,par2);
					Particle[id].Position=Troncature(Particle[id].Position);
					MapDose(false,Dose,Particle[id].Position);
					RendDose(false,Dose,Particle[id].Position);
					alongStepLength+=trueStepLength;
					*totalLength+=trueStepLength;
					bool_cut=false;
				break;
				
				case false:	
					Dose=GlobalLoss(trueStepLength,cutEnergyElectron);	
					GlobalMscScattering(trueStepLength,lengthtoVertex,currentRange,currentLambda,currentEnergy,cutEnergyElectron,par1,par2);
					Particle[id].Position=Troncature(Particle[id].Position);
					MapDose(false,Dose,Particle[id].Position);
					RendDose(false,Dose,Particle[id].Position);
					switch(prc)
					{
						case 0:
							if(bool.verbose==true)
								printf("Ionisation\n");
							GlobalSampleSecondarieElectron(cutEnergyElectron);
							bool_cut=true;											 
						break;
						
						case 1:
							if(bool.verbose==true)
								printf("Bremssthralung\n");
							eSampleSecondarieGamma(cutEnergyGamma);
							bool_cut=true;											 
						break;

						case 2:
							if(bool.verbose==true)
								printf("Annihilation\n");
							Particle[id].Ekin=0.;
						break;	
					}
					alongStepLength=0.;
					*totalLength+=trueStepLength;
				break;						 
			}
			if(bool.ray_tracing==true&&Particle[id].level==0)
				fprintf(Ray_tracing_file,"%lf %lf %lf %d\n",Particle[id].Position.X,Particle[id].Position.Y,Particle[id].Position.Z,Particle[id].nature);
 		}
		if(bool.verbose==true)
			printf("\n");
	}while(Particle[id].Ekin>EkinLimit&&bool_loop==true);
	
	if(Particle[id].Ekin>EkinLimit)
	{
		fragment=GetStepSide(Particle[id].Position,cutStep,&face_cube,true);
		if(bool.verbose==true)
			printf("Fragment %lf\n",fragment);
		currentRange=GetRange(Particle[id].Ekin);	currentEnergy=Particle[id].Ekin;	currentLambda=1./GetLambda(currentEnergy,1);
		trueStepLength=GlobalMscScattering(fragment,cutStep,currentRange,currentLambda,currentEnergy,cutEnergyElectron,par1,par2);
		Particle[id].Position=Troncature(Particle[id].Position);
		*freeLength=alongStepLength+trueStepLength;
		*totalLength+=trueStepLength;
		*face_voxel=face_cube;
		if(bool.ray_tracing==true&&Particle[id].level==0)
				fprintf(Ray_tracing_file,"%lf %lf %lf %d\n",Particle[id].Position.X,Particle[id].Position.Y,Particle[id].Position.Z,Particle[id].nature);
		return	Backscattering(face_cube);
	}
	else
		return	false;
}

void main(int argc, char* argv[])
{
	int	particle,follow_particle,bool_secondaries,bool_ray_tracing,bool_inside,bool_inside_target;
	int	i,j,k,face_voxel,Nbr_part,Dens_part,image,Nbr_images,argument;
	int	indice,idVoxel;
	float	Dose_float;
	double	double_image,double_nbr_images;
	double	height,width,thickness,emission_angle,detection_radius;
	double	Ekin_initiale,EkinLimit,theta,phi;
	double	cutEnergyGamma,cutEnergyElectron;	
	double	DoseTotale,totalLength,freeLength;
	chaine	file_project,file_pTherapy,file_dose,file_dose_bin;

	u64 Time=time(NULL);
	//Time=1405428959;
	randk_seed(Time);
	TimeFile=fopen("Time","w");
	fprintf(TimeFile,"Time %d\n",Time);
	fclose(TimeFile);

	cutEnergyElectron=990.*MeV;	
	//cutEnergyElectron=10.*keV;	
	//cutEnergyElectron=14.0874*keV;	
	//cutEnergyElectron=9.99095*keV;
	//cutEnergyElectron=10.9814*MeV;	
	//cutEnergyElectron=1.*GeV;
	cutEnergyGamma=990.*eV;	
	cutEnergyGamma=10.*GeV;	
	EkinLimit=MinKinEnergy;									 
	xi=1./100.;	
		
	Nbr_part	=	1000000;
	Dens_part	=	00; 
	particle	=	proton;
	Ekin_initiale	=	200.*MeV;
	height			=	200.*mm;
	width				=	200.*mm;
	thickness		=	200.*mm;
	Source.Origine.X		=	0.;	Source.Origine.Y		=	0.;	Source.Origine.Z		=	0.;
	Source.Origine.Z		=	-thickness/2.;
	Source.Direction.X	=	0.;	Source.Direction.Y	=	0.;	Source.Direction.Z	=	1.;
	Detecteur.Spacing		= 1.*mm;	
	Detecteur.Outrange	= 2.*cm;
	emission_angle			=	0./deg;
	detection_radius		=	1.*cm;
	follow_particle			=	particle;
	idsec=MAX_PART-1;
	
	bool.phantom							=	true;
	bool.isotropy							=	false;
	bool.homogenity						=	false;
	bool.source_rotation			=	false;
	bool.hperte								=	true;
	bool.eperte								=	false;
	bool.ebrem								=	false;
	bool.hfluctuation					=	true;
	bool.efluctuation					=	false;
	bool_secondaries					=	false;
	bool.tertiaries						=	false;
	bool.hmultiplescattering	=	true;
	bool.emultiplescattering	=	false;
	bool.lateral							=	true;
	bool.panni								=	false;
	bool.tracking							=	true;
	bool.projection						=	false;
	bool.reconstruction3D			=	false;
	bool.hmapping							=	false;
	bool.emapping							=	false;
	bool.rendement						=	false;
	bool_ray_tracing					=	true;
	bool.verbose							=	false;
	
	if(bool.phantom==true)
	{
		PhantomReader();
		Target.Dimension.X=Voxel.Dimension.X*Voxel.bin.x/2.;
		Target.Dimension.Y=Voxel.Dimension.Y*Voxel.bin.y/2.;
		Target.Dimension.Z=Voxel.Dimension.Z*Voxel.bin.z/2.;
		Volume.Dimension.X=Voxel.Dimension.X/2.;
		Volume.Dimension.Y=Voxel.Dimension.Y/2.;
		Volume.Dimension.Z=Voxel.Dimension.Z/2.;
		Source.Origine.Z=-Target.Dimension.Z;
	}
	else
	{	
		for(i=0;i<MAX_MAT;i++)
			Mater[i]=false;
 		Mater[G4_Water]=true;
    // Mater[Sch_Water]=true;
		// Mater[G4_Bone_Compact]=true;
 		// Mater[G4_Blood]=true;
		Volume.Dimension.X=(width/2.);
		Volume.Dimension.Y=(height/2.);
		Volume.Dimension.Z=(thickness/2.);
		Volume.Center.X=
		Volume.Center.Y=
		Volume.Center.Z=0.;
		Target.Dimension=Volume.Dimension;
		Voxel.Dimension.X=Volume.Dimension.X*2.;
		Voxel.Dimension.X=1.*mm;
		Voxel.Dimension.Y=Volume.Dimension.Y*2.;
		Voxel.Dimension.Y=1.*mm;
		Voxel.Dimension.Z=Volume.Dimension.Z*2.;
		Voxel.Dimension.Z=1.*mm;
		Voxel.bin.x=Voxel.bin.y=Voxel.bin.z=0;
		Voxel.indice.x=Voxel.indice.y=Voxel.indice.z=0;
	}	
	Target.Dimension=Troncature(Target.Dimension);
	Voxel.Dimension=Troncature(Voxel.Dimension);
	Volume.Dimension=Troncature(Volume.Dimension);

	if(bool.hmapping==true||bool.emapping==true)	
	{
		if(bool.phantom==false)
			VoxelValidity();
		MapDose(true,0.,Source.Origine);
	}

	if(bool.rendement==true)
	{
		RendDose(true,0.,Source.Origine);
		Rendement_file=fopen("Rendement.fit","w");
	}	

	Bins					=	1001;
	MinKinEnergy	=	1.*eV;
	MaxKinEnergy	=	250.*MeV;
	if(Bins>MAX_RNG)
		Bins=MAX_RNG;
	
	for(j=0;j<MAX_MAT;j++)
	{
		if(Mater[j]==true)
		{
			id_mat=j;
			Material();
			Generate_table(cutEnergyGamma,cutEnergyElectron);
			// PrintTable(cutEnergyGamma,cutEnergyElectron);
		}
	}
	if(bool.reconstruction3D==true&&argc==1)
		Nbr_images=256;
	else
		Nbr_images=1;
		
	for(image=0;image<Nbr_images;image++)
	{	
		printf("\nProjection %d/%d\n",image+1,Nbr_images);
		double_image			=(double)image;
		double_nbr_images	=(double)Nbr_images;
		if(bool.reconstruction3D==true)
		{
			if(argc>1)
			{
				image=atoi(argv[1]);
				double_image=(double)image;
			}	
			//!§ emission_angle=double_image/double_nbr_images*360./deg;
			emission_angle=double_image/256.*360./deg;
			sprintf(file_project,"../../Fichiers_projection/Project_particle_CPU_%d.fit",image);
			sprintf(file_pTherapy,"../../Fichiers_projection/pTherapy_CPU_%d.fit",image);
		}	
		else
		{
			sprintf(file_project,"Project_particle_CPU.fit");
			sprintf(file_pTherapy,"pTherapy_CPU.fit");
		}

		if(bool.tracking==true)
			Output=fopen(file_pTherapy,"w");
		if(bool.projection==true)
			Projection_file=fopen(file_project,"w");
		if(bool_ray_tracing==true)
			Ray_tracing_file=fopen("../../Python/Particles/Particle.fit","w");
		// FreeOutput=fopen("FreeVar.fit","w");
		DeltaRay=fopen("DeltaRay.fit","w");
		
		if(bool.source_rotation==true)
		{
			if(detection_radius<sqrt(pow(Target.Dimension.X,2.)+pow(Target.Dimension.Z,2.)))
				detection_radius=sqrt(pow(Target.Dimension.X,2.)+pow(Target.Dimension.Z,2.));
			Detecteur.Dimension.Z=detection_radius;
			Source.Origine.Z=-detection_radius;
		}
		else
		{	
			//Detecteur.Dimension.Z=sqrt(pow(Target.Dimension.X,2.)+pow(Target.Dimension.Z,2.));;
			Detecteur.Dimension.Z=Target.Dimension.Z+Detecteur.Spacing;
			//Source.Origine.Z=-Detecteur.Dimension.Z;
			emission_angle=0.*deg;	
		}
		Detecteur.Dimension.X=Target.Dimension.X*fabs(cos(emission_angle))+Target.Dimension.Z*fabs(sin(emission_angle))+Detecteur.Outrange;
		Detecteur.Dimension.Y=Target.Dimension.Y+Detecteur.Outrange;	
		if(bool.homogenity==true)
		{	
			Source.Dimension.X=Detecteur.Dimension.X-Detecteur.Outrange;
			Source.Dimension.Y=Detecteur.Dimension.Y-Detecteur.Outrange;
		}
		else
		{
			Source.Dimension.X=0.*cm;
			Source.Dimension.Y=0.*cm;
		}	
		Source.Dimension.Z=Source.Origine.Z;
		Source.Angle=emission_angle;
		if(bool.homogenity==true&&Dens_part>0)
			Nbr_part=Dens_part*4.*Source.Dimension.X*Source.Dimension.Y;	
		
		if(bool_ray_tracing==true)
			PrintGeometry();
				
		// id_mat=Sch_Water;
		id_mat=G4_Water;
		// id_mat=G4_Bone_Compact;

		for(i=0;i<Nbr_part;i++)							 	 
		{
			Time++;
			randk_seed(Time);

			if(i==0)
			{	
				printf("Première particule des %d évenements\n",Nbr_part);
				bool.ray_tracing=bool_ray_tracing;
			}
			else
			{
				bool.ray_tracing=false;
				bool.verbose=false;
			}

			Initialize(particle,Ekin_initiale);
			bool.secondaries=bool_secondaries;
			if(bool.secondaries==false)
				bool.tertiaries=false;
			for(id=0;id<=idsec;id++)
			{
				Particle[id].Position=Troncature(Particle[id].Position);
				// if(id==1)	
				// 	bool.verbose=true;
				// else
				// 	bool.verbose=false;

				if(bool.verbose==true)
				{
					if(id==0)
						printf("\n#############		Début du primaire		#############\n\n");
					printf("Début pour la particule : %d du run : %d; Energie : %lf; Level : %d\n",id,i,Particle[id].Ekin,Particle[id].level);
					printf("Position : ");
					PrintVector(Particle[id].Position);
					printf("Moment : ");
					PrintVector(Particle[id].Momentum);
				}
				if(bool.phantom==true)
					FirstVoxel(Particle[id].Position,Time);
				Parameters();
				totalLength=0.*mm;
				freeLength=0.*mm;
				do{	
					do{
						bool_inside=Navigate(&totalLength,&freeLength,cutEnergyElectron,cutEnergyGamma,EkinLimit,&face_voxel);	
					}while(bool_inside==true);
					if(bool.verbose==true)
						printf("Sortie du voxel\n");
					if(Particle[id].Ekin>EkinLimit)
						bool_inside_target=NextVoxel(face_voxel);
				}while(Particle[id].Ekin>EkinLimit&&bool_inside_target==true);	
				if(bool.verbose==true&&id==0)
					printf("\n#############		Fin du primaire		#############\n\n");
				if(Particle[id].Ekin>EkinLimit&&Particle[id].nature==follow_particle)
				{		
					Particle[id].Position=Troncature(Particle[id].Position);
					if(DetectorTrajectory()==true)
					{
						//theta=Theta();	phi=Phi();
						if(bool.source_rotation==true)
							AngularCorrection();
						if(bool.projection==true)
							// fprintf(Projection_file,"%8.6lf	%8.6lf	%8.6lf	%8.7lf	%8.7lf	%8.7lf	%8.7lf	%8.6lf	%8.6lf	%8.6lf	%8.7lf	%8.7lf	%8.7lf\n",
							// 				Source.Origine.X/mm,Source.Origine.Y/mm,Source.Origine.Z/mm,Source.Direction.X,Source.Direction.Y,Source.Direction.Z,Particle[id].Ekin,
							// 				Particle[id].Position.X/mm,Particle[id].Position.Y/mm,totalLength/mm,Particle[id].Momentum.X,Particle[id].Momentum.Y,Particle[id].Momentum.Z);	
							fprintf(Projection_file,"%5.5lf	%5.5lf	%5.5lf	%5.5lf	%5.5lf	%5.5lf	%5.5lf	%5.5lf\n",
											Source.Origine.X/mm,Source.Origine.Y/mm,Particle[id].Ekin,
											Particle[id].Position.X/mm,Particle[id].Position.Y/mm,
											Particle[id].Momentum.X,Particle[id].Momentum.Y,Particle[id].Momentum.Z);	
						if(bool.tracking==true)
							fprintf(Output,"%8.7lf	%8.6lf	%8.6lf	%8.6lf	%8.7lf	%8.7lf	%8.7lf\n",
											Particle[id].Ekin,
											Particle[id].Position.X/mm,Particle[id].Position.Y/mm,Particle[id].Position.Z/mm,//totalLength/mm,
											Particle[id].Momentum.X,Particle[id].Momentum.Y,Particle[id].Momentum.Z);
							//fprintf(Output,"%8.7lf\n",Particle[id].Ekin);
						Particle[id].Ekin=0.;
					}
				}
				// idsec=0; 
			}	
			if((i+1)%500000==0)
				printf("Evenement	%d/%d\n",i+1,Nbr_part);
		}
		if(bool.hmapping==true||bool.emapping==true)
		{
			if(argc>1)
			{
				argument=atoi(argv[1]);
				sprintf(file_dose,"../../Fichiers_dose/Totale_dosimetry_CPU_%d.fit",argument);
				sprintf(file_dose_bin,"../../Fichiers_dose/Totale_dosimetry_CPU_%d.bin",argument);
			}	
			else
			{
				sprintf(file_dose,"Totale_dosimetry_CPU.fit");
				sprintf(file_dose_bin,"Totale_dosimetry_CPU.bin");
			}
			Dosimetry=fopen(file_dose,"w");
			Dosimetry_bin=fopen(file_dose_bin,"wb");
			DoseTotale=0.;
			for(k=0;k<Voxel.bin.z;k++)	for(j=0;j<Voxel.bin.y;j++)	for(i=0;i<Voxel.bin.x;i++)	
			{
				idVoxel=k*Voxel.bin.x*Voxel.bin.y+j*Voxel.bin.x+i;
				if(Vecteur_Dose_map[idVoxel]>0.)
				{
					fprintf(Dosimetry,"%d	%d	%d	%E\n",i,j,k,Vecteur_Dose_map[idVoxel]);
					DoseTotale+=Vecteur_Dose_map[idVoxel];
				}
			}
			for(k=0;k<Voxel.bin.z;k++)	for(j=0;j<Voxel.bin.y;j++) for(i=0;i<Voxel.bin.x;i++)	
			{
				idVoxel=k*Voxel.bin.x*Voxel.bin.y+j*Voxel.bin.x+i;
				Dose_float=(float)Vecteur_Dose_map[idVoxel];
				fwrite(&Dose_float,sizeof(float),1,Dosimetry_bin);
			}			
			printf("Energie totale déposée %5.2lf GeV soit %5.2lf MeV/part\n",EnergyLost/GeV,EnergyLost/Nbr_part);
			printf("Dose totale déposée %5.2E Gy soit %5.2E Gy/part\n",DoseTotale,DoseTotale/Nbr_part);
			fclose(Dosimetry);
			fclose(Dosimetry_bin);
			free(Vecteur_Dose_map);
		}				
		if(bool.rendement==true)
		{
			for(i=0;i<MAX_RNG;i++)
				if(RendementMapping[i]>0.)
					fprintf(Rendement_file,"%d %lf\n",i,RendementMapping[i]);
			fclose(Rendement_file);
		}	
		if(bool.phantom==true)
			free(Vecteur_Voxel_mat);
		if(bool.tracking==true)
			fclose(Output);
		if(bool.projection==true)
			fclose(Projection_file);
		if(bool_ray_tracing==true)
			fclose(Ray_tracing_file);	
		// fclose(FreeOutput);	
		fclose(DeltaRay);
	}	
}