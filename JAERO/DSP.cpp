//---------------------------------------------------------------------------


#include "DSP.h"
#include <QDebug>

//---------------------------------------------------------------------------


float arctan2_fast_maybe(cpx_type point)
{
   double y=point.imag();
   double x=point.real();
   double angle;
   double abs_y = fabs(y)+1e-10;
   if (x>=0)
   {
      double r = (x - abs_y) / (x + abs_y);
      angle=0.1963 * powf(r,3.0) - 0.9817 * r + M_PI/4.0;
   }
    else
    {
      double r = (x + abs_y) / (abs_y - x);
      angle=0.1963 * powf(r,3.0) - 0.9817 * r + 3.0*M_PI/4.0;
    }
   if (y < 0) return(-angle);
    else return(angle);
}

//this is a problem. If WaveTableComplex is used like WaveTableComplex wt; in global then strange things happen like needing to init wt_cis twice. This doesn't happen if used in a class or as a pointer and created later.
std::vector<cpx_type> WaveTableComplex::wt_cis;

TrigLookUp tringlookup;

TrigLookUp::TrigLookUp()
{
    //--load wavetables (trig)
    SinWT.resize(WTSIZE);
    CosWT.resize(WTSIZE);
    CosINV.resize(WTSIZE);
    CISWT.resize(WTSIZE);
    int cosine,i;
    for(i=0;i<WTSIZE;i++){SinWT[i]=(sin(2*M_PI*((double)i)/WTSIZE));}
    for(i=0;i<WTSIZE;i++){CosWT[i]=(sin(M_PI_2+2*M_PI*((double)i)/WTSIZE));}
    for(i=0;i<WTSIZE;i++)
    {
        CISWT[i]=cpx_type(CosWT[i],SinWT[i]);
        cosine=(2.0*((double)i)/(WTSIZE-1.0))-1.0;
        if (cosine>1)cosine=1;
        else if(cosine<-1)cosine=-1;
        CosINV[i]=acos(cosine);
    }
    //--
}

WaveTable::WaveTable()
{
    last_WTptr=0;
    samplerate=48000;
    freq=1000;
    WTstep=(1000.0)*WTSIZE/(48000);//default
    WTptr=0;
    FractionOfSampleItPassesBy=0.0;
}

WaveTable::WaveTable(int _freq,int _samplerate)
{
    last_WTptr=0;
    freq=_freq;
    samplerate=_samplerate;
    if(freq<0)freq=0;
    WTstep=((double)freq)*WTSIZE/((float)samplerate);
    WTptr=0;
    FractionOfSampleItPassesBy=0.0;
}

WaveTable::~WaveTable()
{

}

double  WaveTable::DistancebetweenWT(double WTptr1, double WTptr2)
{
    double dif1=WTptr1-WTptr2;
    double dif2=WTptr2-WTptr1;
    if(dif1<0)dif1+=WTSIZE;
    if(dif2<0)dif2+=WTSIZE;
    double dif=dif1;
    if(dif2<dif1)dif=-dif2;
    dif=dif/WTSIZE;
    return dif;
}

void  WaveTable::WTnextFrame()
{
    JASSERT(WTptr>=0.0);
    if(WTstep<0)WTstep=0;
    last_WTptr=WTptr;
    WTptr+=WTstep;
    while(((int)WTptr)>=WTSIZE) WTptr-=WTSIZE;
}

cpx_type  WaveTable::WTCISValue()
{
    tint=(int)WTptr;
    if(tint>=WTSIZE)tint=0;
    if(tint<0)tint=WTSIZE-1;
    return tringlookup.CISWT[tint];
}

double   WaveTable::WTSinValue()
{
    tint=(int)WTptr;
    if(tint>=WTSIZE)tint=0;
    if(tint<0)tint=WTSIZE-1;
    return tringlookup.SinWT[tint];
}

double   WaveTable::WTSinValue(double PlusFractOfSample)
{
    double ts=WTptr;
    ts+=PlusFractOfSample*WTstep;
    if(((int)ts)>=WTSIZE) ts-=WTSIZE;
    tint=(int)ts;
    if(tint>=WTSIZE)tint=0;
    if(tint<0)tint=WTSIZE-1;
    return tringlookup.SinWT[tint];
}

double  WaveTable::WTCosValue()
{
    tint=(int)WTptr;
    if(tint>=WTSIZE)tint=0;
    if(tint<0)tint=WTSIZE-1;
    return tringlookup.CosWT[tint];
}

double   WaveTable::WTCosValue(double PlusFractOfSample)
{
    double ts=WTptr;
    ts+=PlusFractOfSample*WTstep;
    if(((int)ts)>=WTSIZE) ts-=WTSIZE;
    tint=(int)ts;
    if(tint>=WTSIZE)tint=0;
    if(tint<0)tint=WTSIZE-1;
    return tringlookup.CosWT[tint];
}

void  WaveTable::WTsetFreq(int _freq,int _samplerate)
{
    samplerate=_samplerate;
    freq=_freq;
    if(freq<0)freq=0;
    WTstep=((double)freq)*WTSIZE/((float)samplerate);
    while(((int)WTptr)>=WTSIZE) WTptr-=WTSIZE;
}

void  WaveTable::SetFreq(double _freq,int _samplerate)
{
    freq=_freq;
    samplerate=_samplerate;
    if(freq<0)freq=0;
    WTstep=(freq)*((double)WTSIZE)/((float)samplerate);
    while(((int)WTptr)>=WTSIZE) WTptr-=WTSIZE;
}

void  WaveTable::SetFreq(double _freq)
{
    freq=_freq;
    if(freq<0)freq=0;
    WTstep=(freq)*((double)WTSIZE)/samplerate;
}

double WaveTable::GetFreqHz()
{
    return freq;
}

void  WaveTable::IncreseFreqHz(double freq_hz)
{
    freq_hz+=freq;
    SetFreq(freq_hz);
}

void  WaveTable::IncresePhaseDeg(double phase_deg)
{
    phase_deg+=(360.0*WTptr/((double)WTSIZE));
    SetPhaseDeg(phase_deg);
}

void  WaveTable::SetPhaseDeg(double phase_deg)
{
    phase_deg=std::fmod(phase_deg,360.0);
    while(phase_deg<0)phase_deg+=360.0;
    WTptr=(phase_deg/360.0)*((double)WTSIZE);
}

void  WaveTable::SetPhaseCycles(double phase_cycles)//untested
{
    phase_cycles=std::fmod(phase_cycles,1.0);
    while(phase_cycles<0)phase_cycles+=1.0;
    assert(phase_cycles>=0);
    assert(phase_cycles<1);
    WTptr=phase_cycles*((double)WTSIZE);
}


double WaveTable::GetPhaseDeg()
{
    return (360.0*WTptr/((double)WTSIZE));
}

double WaveTable::GetFreqTest()
{
    return WTstep;
}

bool  WaveTable::IfPassesPointNextTime()
{
    FractionOfSampleItPassesBy=WTptr+(WTstep-WTSIZE);
    if(FractionOfSampleItPassesBy<0)return false;
    FractionOfSampleItPassesBy/=WTstep;
    return true;
}

bool  WaveTable::IfPassesPointNextTime(double FractionOfWave)
{
    FractionOfSampleItPassesBy=(FractionOfWave*WTSIZE)-WTptr;
    if(FractionOfSampleItPassesBy<0)FractionOfSampleItPassesBy+=WTSIZE;
    if(FractionOfSampleItPassesBy<WTstep)
    {
        FractionOfSampleItPassesBy=(WTstep-FractionOfSampleItPassesBy)/WTstep;
        return true;
    }
    return false;
}

bool  WaveTable::IfHavePassedPoint(double FractionOfWave)
{
    double t_last_WTptr=last_WTptr;
    double t_WTptr=WTptr;
    double pt=(FractionOfWave*WTSIZE);
    t_last_WTptr-=pt;
    t_WTptr-=pt;
    if(t_last_WTptr<0.0)t_last_WTptr+=WTSIZE;
    if(t_WTptr<0.0)t_WTptr+=WTSIZE;
    if((t_last_WTptr>WTSIZE_3_4)&&(t_WTptr<WTSIZE_1_4))
    {
        //FractionOfSampleItPassesBy=(WTstep-t_WTptr)/WTstep;
        FractionOfSampleItPassesBy=t_WTptr/WTstep;//WTSIZE;
        return true;
    }
    return false;
}


bool  WaveTable::IfPassesPointNextTime_frames(double NumberOfFrames)
{
    FractionOfSampleItPassesBy=NumberOfFrames-WTptr;
    if(FractionOfSampleItPassesBy<0)FractionOfSampleItPassesBy+=WTSIZE;
    if(FractionOfSampleItPassesBy<WTstep)
    {
        FractionOfSampleItPassesBy=(WTstep-FractionOfSampleItPassesBy)/WTstep;
        return true;
    }
    return false;
}





void  WaveTable::SetWTptr(double FractionOfWave,double PlusNumberOfFrames)
{
    while(FractionOfWave>=1)FractionOfWave-=1;
    while(FractionOfWave<0)FractionOfWave+=1;
    WTptr=(FractionOfWave*WTSIZE);
    WTptr+=PlusNumberOfFrames*WTstep;
    while(((int)WTptr)>=WTSIZE)WTptr-=WTSIZE;
    while(((int)WTptr)<0)WTptr+=WTSIZE;
}

//---------------------------------------------------------------------------



FIR::FIR(int _NumberOfPoints)
{
        int i;
        points=0;buff=0;

        NumberOfPoints=_NumberOfPoints;
        buffsize=NumberOfPoints+1;
        points=new double[NumberOfPoints];
        for(i=0;i<NumberOfPoints;i++)points[i]=0;
        buff=new double[buffsize];
        for(i=0;i<buffsize;i++)buff[i]=0;
        ptr=0;
        outsum=0;
}

FIR::~FIR()
{
        if(points)delete [] points;
        if(buff)delete [] buff;
}

double  FIR::FIRUpdateAndProcess(double sig)
{
        buff[ptr]=sig;
        ptr++;if(ptr>=buffsize)ptr=0;//ptr%=buffsize;
        int tptr=ptr;
        outsum=0;
        for(int i=0;i<NumberOfPoints;i++)
        {
                outsum+=points[i]*buff[tptr];
                tptr++;if(tptr>=buffsize)tptr=0;//tptr%=buffsize;
        }
        return outsum;
}

void  FIR::FIRUpdate(double sig)
{
        buff[ptr]=sig;
        ptr++;ptr%=buffsize;
}

double  FIR::FIRProcess(double FractionOfSampleOffset)
{
        double nextp=FractionOfSampleOffset;
        double thisp=1-nextp;
        int tptr=ptr;
        int nptr=ptr;
        outsum=0;
        nptr++;nptr%=buffsize;
        for(int i=0;i<NumberOfPoints;i++)
        {
                outsum+=points[i]*(buff[tptr]*thisp+buff[nptr]*nextp);
                tptr=nptr;
                nptr++;nptr%=buffsize;
        }
        return outsum;
}

double  FIR::FIRUpdateAndProcess(double sig, double FractionOfSampleOffset)
{
        buff[ptr]=sig;
        ptr++;ptr%=buffsize;
        double nextp=FractionOfSampleOffset;
        double thisp=1-nextp;
        int tptr=ptr;
        int nptr=ptr;
        outsum=0;
        nptr++;nptr%=buffsize;
        for(int i=0;i<NumberOfPoints;i++)
        {
                outsum+=points[i]*(buff[tptr]*thisp+buff[nptr]*nextp);
                tptr=nptr;
                nptr++;nptr%=buffsize;
        }
        return outsum;
}

void  FIR::FIRSetPoint(int point, double value)
{
    JASSERT(point>=0);
    JASSERT(point<NumberOfPoints);
    if((point<0)||(point>=NumberOfPoints))return;
    points[point]=value;
}

//-----------------
AGC::AGC(double _SecondsToAveOver,double _Fs)
{
    JASSERT(_Fs>1);
    JASSERT(_Fs<1000000);
    AGCMASz=round(_SecondsToAveOver*_Fs);
    JASSERT(AGCMASz>0);
    AGCMASum=0;
    AGCMABuffer=new double[AGCMASz];
    for(int i=0;i<AGCMASz;i++)AGCMABuffer[i]=0;
    AGCMAPtr=0;
    AGCVal=0;
}

double AGC::Update(double sig)
{
    AGCMASum=AGCMASum-AGCMABuffer[AGCMAPtr];
    AGCMASum=AGCMASum+fabs(sig);
    AGCMABuffer[AGCMAPtr]=fabs(sig);
    AGCMAPtr++;AGCMAPtr%=AGCMASz;
    AGCVal=1.414213562/fmax(AGCMASum/((double)AGCMASz),0.000001);
    AGCVal=fmax(AGCVal,0.000001);
    return AGCVal;
}

AGC::~AGC()
{
    if(AGCMASz)delete [] AGCMABuffer;
}

//---------------------

MovingAverage::MovingAverage(int number)
{
    //MASz=round(number);
    MASz=number;
    JASSERT(MASz>0);
    MASum=0;
    MABuffer=new double[MASz];
    for(int i=0;i<MASz;i++)MABuffer[i]=0;
    MAPtr=0;
    Val=0;
}

void MovingAverage::Zero()
{
    for(int i=0;i<MASz;i++)MABuffer[i]=0;
    MAPtr=0;
    Val=0;
    MASum=0;
}

double MovingAverage::Update(double sig)
{
    MASum=MASum-MABuffer[MAPtr];
    MASum=MASum+fabs(sig);
    MABuffer[MAPtr]=fabs(sig);
    MAPtr++;MAPtr%=MASz;
    Val=MASum/((double)MASz);
    return Val;
}

double MovingAverage::UpdateSigned(double sig)
{
    MASum=MASum-MABuffer[MAPtr];
    MASum=MASum+(sig);
    MABuffer[MAPtr]=(sig);
    MAPtr++;MAPtr%=MASz;
    Val=MASum/((double)MASz);
    return Val;
}

MovingAverage::~MovingAverage()
{
    if(MASz)delete [] MABuffer;
}
//---------------------

MSEcalc::MSEcalc(int number)
{
    pointmean = new MovingAverage(number);
    msema = new MovingAverage(number);
    mse=0;
}
MSEcalc::~MSEcalc()
{
    delete pointmean;
    delete msema;
}
void MSEcalc::Zero()
{
    pointmean->Zero();
    msema->Zero();
    mse=0;
}
double MSEcalc::Update(cpx_type pt_qpsk)
{
    double tda,tdb;
    cpx_type tcpx;
    pointmean->Update(std::abs(pt_qpsk));
    double mu=pointmean->Val;
    if(mu<0.000001)mu=0.000001;
    tcpx=sqrt(2)*pt_qpsk/mu;
    tda=(fabs(tcpx.real())-1.0);
    tdb=(fabs(tcpx.imag())-1.0);
    mse=msema->Update((tda*tda)+(tdb*tdb));
    return mse;
}

//---------------------
MovingVar::MovingVar(int number)
{
    E = new MovingAverage(number);
    E2 = new MovingAverage(number);
}

double MovingVar::Update(double sig)
{
    E2->Update(sig*sig);
    Mean=E->Update(sig);
    Var=(E2->Val)-(E->Val*E->Val);
    return Var;
}

MovingVar::~MovingVar()
{
    delete E;
    delete E2;
}
//---------------------

MSKEbNoMeasure::MSKEbNoMeasure(int number)
{
    E = new MovingAverage(number);
    E2 = new MovingAverage(number);
}

double MSKEbNoMeasure::Update(double sig)
{
    E2->Update(sig*sig);
    Mean=E->Update(sig);
    Var=(E2->Val)-(E->Val*E->Val);
    double alpha=sqrt(2)/Mean;
    //double tebno=10.0*(log10(2.0)-log10(((Var*alpha*alpha)- 0.0085 ))+log10(1.5));//0.0085 for the non constant modulus after the matched filter
    double tebno=10.0*(log10(2.0)-log10(((Var*alpha*alpha)- 0.0085 )))-5.0;// this one matches matlab better
    if(std::isnan(tebno))tebno=50;
    if(tebno>50.0)tebno=50;
    EbNo=EbNo*0.8+0.2*tebno;
    return EbNo;
}

MSKEbNoMeasure::~MSKEbNoMeasure()
{
    delete E;
    delete E2;
}


//----------------------


DiffDecode::DiffDecode()
{
    laststate=false;
    lastsoftstate= -1;
}

bool DiffDecode::Update(bool state)
{
    bool res=(state+laststate)%2;
    laststate=state;
    return res;
}


double DiffDecode::UpdateSoft(double soft)
{
    double retval = 0;

    // if the previous value is a zero and the current also zero just return zero
    if(soft < 0 && lastsoftstate < 0)
    {

        // last value is negative so just return to indicate a zero
        retval = lastsoftstate;
        lastsoftstate = soft;
    }

     // if the previous value is one and the current also one
     else if(soft > 0 && lastsoftstate > 0)
     {

        // last value is postive so flip sign to indicate zero
        retval =- lastsoftstate;
        lastsoftstate = soft;
     }
      else
      {


        // retval and soft have different signs, so always return positive
        retval = std::fabs(lastsoftstate);
        lastsoftstate = soft;
      }

    return retval;

}

//-----------------------

bool  BaceConverter::LoadSymbol(int Value)
{
        if(Value<0)
        {
                Value=0;
                ErasurePtr=BuffPosPtr+NumberOfBitsForInSize;
        }

        Buff|=(Value<<BuffPosPtr);
        BuffPosPtr+=NumberOfBitsForInSize;
        if(BuffPosPtr>=NumberOfBitsForOutSize)DataAvailable=true;
         else DataAvailable=false;
        return DataAvailable;
}

bool  BaceConverter::GetNextSymbol()
{
         if(BuffPosPtr>=NumberOfBitsForOutSize)
         {
                if(ErasurePtr>0)
                {
                        Result=-1;
                        ErasurePtr-=NumberOfBitsForOutSize;
                }
                 else Result=Buff&OutMaxVal;
                Buff=(Buff>>NumberOfBitsForOutSize);
                BuffPosPtr-=NumberOfBitsForOutSize;
         }
         if(BuffPosPtr>=NumberOfBitsForOutSize)DataAvailable=true;
          else DataAvailable=false;
         return DataAvailable;
}

 BaceConverter::BaceConverter()
{
        SetInNumberOfBits(4);
        SetOutNumberOfBits(8);
}


void  BaceConverter::SetInNumberOfBits(int NumberOfBits)
{
        InMaxVal=((int)pow(2,NumberOfBits))-1;
        NumberOfBitsForInSize=NumberOfBits;
        Reset();
}

void  BaceConverter::SetOutNumberOfBits(int NumberOfBits)
{
        OutMaxVal=((int)pow(2,NumberOfBits))-1;
        NumberOfBitsForOutSize=NumberOfBits;
        Reset();
}


void  BaceConverter::Reset()
{
        Result=-1;
        DataAvailable=false;
        Buff=0;
        BuffPosPtr=0;
        Result=0;
        ErasurePtr=0;
}

//-----

IIR::IIR()
{
    a.resize(3);
    b.resize(3);
    b[0]=0.00032714218939589035;
    b[1]=0;
    b[2]=0.00032714218939589035;
    a[0]=1;
    a[1]=-0.39005299948210803;
    a[2]= 0.99934571562120822;
    init();
}

void IIR::init()
{
    JASSERT(a.size()>=1);
    JASSERT(b.size()>=1);
    buff_x.resize(b.size());
    buff_y.resize(a.size()-1);
    buff_x_ptr=0;
    buff_y_ptr=0;
    buff_x.fill(0);
    buff_y.fill(0);
}

double  IIR::update(double sig)
{

    buff_x.resize(b.size());
    buff_y.resize(a.size()-1);
    if(buff_x_ptr>=buff_x.size())buff_x_ptr=0;
    if(buff_y_ptr>=buff_y.size())buff_y_ptr=0;

    ASSERTCH(buff_x,buff_x_ptr);
    buff_x[buff_x_ptr]=sig;
    buff_x_ptr++;buff_x_ptr%=buff_x.size();

    y=0;

    //int tp=buff_x_ptr;

    for(int i=b.size()-1;i>=0;i--)
    {
        ASSERTCH(buff_x,buff_x_ptr);
        ASSERTCH(b,i);
        y+=buff_x[buff_x_ptr]*b[i];
        buff_x_ptr++;buff_x_ptr%=buff_x.size();
    }

    for(int i=a.size()-1;i>=1;i--)
    {
        ASSERTCH(buff_y,buff_y_ptr);
        ASSERTCH(a,i);
        y-=buff_y[buff_y_ptr]*a[i];
        buff_y_ptr++;buff_y_ptr%=buff_y.size();
    }

    ASSERTCH(a,0);
    y/=a[0];

    ASSERTCH(buff_y,buff_y_ptr);
    buff_y[buff_y_ptr]=y;
    buff_y_ptr++;buff_y_ptr%=buff_y.size();

    return y;


/*            a(1)*y(n) = b(1)*x(n) + b(2)*x(n-1) + ... + b(nb+1)*x(n-nb)
                - a(2)*y(n-1) - ... - a(na+1)*y(n-na)*/
}

//----------------

//for alpha==1
//not sure if correct for any other Fs other than 48000 and fb 10500
OQPSKEbNoMeasure::OQPSKEbNoMeasure(int number,double _Fs,double _fb)
{
    Fs=_Fs;
    fb=_fb;
    E = new MovingAverage(number);
    E2 = new MovingAverage(number);
}

double OQPSKEbNoMeasure::Update(double sig)
{
    E2->Update(sig*sig);
    Mean=E->Update(sig);
    double MeanSquared=Mean*Mean;
    Var=(E2->Val)-(E->Val*E->Val);
    Var-=(0.024709*MeanSquared);//remove non constant of OQPSK after RRC
    double mvr=(((Fs*MeanSquared/(2.0*fb*Var)))*0.13743);//calibrated using matlab
    if(mvr<0.000000001)mvr=0.000000001;
    double tebno=10.0*log10(mvr);
    if(std::isnan(tebno))tebno=50;
    if(tebno>50.0)tebno=50;
    if(tebno<0.0)tebno=0;//no real reason but there is no way that less will be usful
    EbNo=EbNo*0.8+0.2*tebno;
    return EbNo;
}

OQPSKEbNoMeasure::~OQPSKEbNoMeasure()
{
    delete E;
    delete E2;
}



//---fast FIR (real complex and hilbert)

double sinc_normalized(double val)
{
    if (val==0)return 1.0;
    return (sin(M_PI*val)/(M_PI*val));
}

QVector<kiss_fft_scalar> QJFilterDesign::LowPassHanning(double FrequencyCutOff, double SampleRate, int Length)
{
    QVector<kiss_fft_scalar> h;
    if(Length<1)return h;
    if(!(Length%2))Length++;
    int j=1;
    for(int i=(-(Length-1)/2);i<=((Length-1)/2);i++)
    {
        double w=0.5*(1.0-cos(2.0*M_PI*((double)j)/((double)(Length))));
        kiss_fft_scalar x;
        x=(w*(2.0*FrequencyCutOff/SampleRate)*sinc_normalized(2.0*FrequencyCutOff*((double)i)/SampleRate));
        h.push_back(x);
        j++;
    }

    return h;

/* in matlab this function is
idx = (-(Length-1)/2:(Length-1)/2);
hideal = (2*FrequencyCutOff/SampleRate)*sinc(2*FrequencyCutOff*idx/SampleRate);
h = hanning(Length)' .* hideal;
*/

}

QVector<kiss_fft_scalar> QJFilterDesign::HighPassHanning(double FrequencyCutOff, double SampleRate, int Length)
{
    QVector<kiss_fft_scalar> h;
    if(Length<1)return h;
    if(!(Length%2))Length++;

    QVector<kiss_fft_scalar> h1;
    QVector<kiss_fft_scalar> h2;
    h2.fill(0,Length);
    h2[(Length-1)/2]=1;
    h1=LowPassHanning(FrequencyCutOff,SampleRate,Length);
    if((h1.size()==Length)&&(h2.size()==Length))
    {
        kiss_fft_scalar val;
        for(int i=0;i<Length;i++)
        {
            val=h2[i]-h1[i];
            h.push_back(val);
        }
    }

    return h;
}

QVector<kiss_fft_scalar> QJFilterDesign::BandPassHanning(double LowFrequencyCutOff,double HighFrequencyCutOff, double SampleRate, int Length)
{
    QVector<kiss_fft_scalar> h;
    if(Length<1)return h;
    if(!(Length%2))Length++;

    QVector<kiss_fft_scalar> h1;
    QVector<kiss_fft_scalar> h2;

    h2=LowPassHanning(HighFrequencyCutOff,SampleRate,Length);
    h1=LowPassHanning(LowFrequencyCutOff,SampleRate,Length);

    if((h1.size()==Length)&&(h2.size()==Length))
    {
        kiss_fft_scalar val;
        for(int i=0;i<Length;i++)
        {
            val=h2[i]-h1[i];
            h.push_back(val);
        }
    }

    return h;
}

QJFastFIRFilter::QJFastFIRFilter(QObject *parent) : QObject(parent)
{
    QVector<kiss_fft_cpx> tvect;
    kernelsize=0;
    kiss_fft_cpx asample;
    asample.r=1;
    asample.i=0;
    tvect.push_back(asample);
    nfft=2;
    cfg=kiss_fastfir_alloc_complex(tvect.data(),tvect.size(),&nfft,0,0);
    reset();
}

int QJFastFIRFilter::setKernel(QVector<kiss_fft_scalar> imp_responce)
{
    int _nfft=imp_responce.size()*4;//rule of thumb
    _nfft=pow(2.0,(ceil(log2(_nfft))));
    return setKernel(imp_responce,_nfft);
}

int QJFastFIRFilter::setKernel(QVector<kiss_fft_cpx> imp_responce)
{
    int _nfft=imp_responce.size()*4;//rule of thumb
    _nfft=pow(2.0,(ceil(log2(_nfft))));
    return setKernel(imp_responce,_nfft);
}

int QJFastFIRFilter::setKernel(QVector<cpx_type> imp_responce)
{
    int _nfft=imp_responce.size()*4;//rule of thumb
    _nfft=pow(2.0,(ceil(log2(_nfft))));
    return setKernel(imp_responce,_nfft);
}

int QJFastFIRFilter::setKernel(QVector<kiss_fft_cpx> imp_responce,int _nfft)
{
    if(!imp_responce.size())return nfft;
    free(cfg);
    _nfft=pow(2.0,(ceil(log2(_nfft))));
    nfft=_nfft;
    cfg=kiss_fastfir_alloc_complex(imp_responce.data(),imp_responce.size(),&nfft,0,0);
    reset();
    kernelsize=imp_responce.size();
    return nfft;
}

int QJFastFIRFilter::setKernel(QVector<cpx_type> imp_responce,int _nfft)
{
    if(!imp_responce.size())return nfft;
    free(cfg);
    _nfft=pow(2.0,(ceil(log2(_nfft))));
    nfft=_nfft;
    QVector<kiss_fft_cpx> timp_responce;
    timp_responce.resize(imp_responce.size());
    for(int i=0;i<imp_responce.size();i++)
    {
        timp_responce[i].r=imp_responce[i].real();
        timp_responce[i].i=imp_responce[i].imag();
    }
    cfg=kiss_fastfir_alloc_complex(timp_responce.data(),timp_responce.size(),&nfft,0,0);
    reset();
    kernelsize=timp_responce.size();
    return nfft;
}

int QJFastFIRFilter::setKernel(QVector<kiss_fft_scalar> imp_responce,int _nfft)
{
    if(!imp_responce.size())return nfft;
    free(cfg);
    _nfft=pow(2.0,(ceil(log2(_nfft))));
    nfft=_nfft;
    QVector<kiss_fft_cpx> timp_responce;
    timp_responce.resize(imp_responce.size());
    for(int i=0;i<imp_responce.size();i++)
    {
        timp_responce[i].r=imp_responce[i];
        timp_responce[i].i=0;
    }
    cfg=kiss_fastfir_alloc_complex(timp_responce.data(),timp_responce.size(),&nfft,0,0);
    reset();
    kernelsize=timp_responce.size();
    return nfft;
}

int QJFastFIRFilter::updateKernel(const QVector<kiss_fft_cpx> &imp_responce)
{
    assert(imp_responce.size()>0);
    if(kernelsize!=imp_responce.size())
    {
        return setKernel(imp_responce);
    }

    size_t nfft=cfg->nfft;
    size_t  n_imp_resp=imp_responce.size();
    const kiss_fft_cpx * imp_resp=imp_responce.data();

    assert((((int)nfft) - ((int)n_imp_resp) + 1)>=0);

    /*zero pad in the middle to left-rotate the impulse response
      This puts the scrap samples at the end of the inverse fft'd buffer */
    cfg->tmpbuf[0] = imp_resp[ n_imp_resp - 1 ];
    for (size_t i=0;i<n_imp_resp - 1; ++i) {
        cfg->tmpbuf[ nfft - n_imp_resp + 1 + i ] = imp_resp[ i ];
    }

    kiss_fft(cfg->fftcfg,cfg->tmpbuf,cfg->fir_freq_resp);

    /* TODO: this won't work for fixed point */
    float scale = 1.0 / ((float)cfg->nfft);

    for (size_t i=0; i < cfg->n_freq_bins; ++i ) {
    #ifdef USE_SIMD
        cfg->fir_freq_resp[i].r *= _mm_set1_ps(scale);
        cfg->fir_freq_resp[i].i *= _mm_set1_ps(scale);
    #else
        cfg->fir_freq_resp[i].r *= scale;
        cfg->fir_freq_resp[i].i *= scale;
    #endif
    }


    return nfft;


}

int QJFastFIRFilter::updateKernel(const QVector<cpx_type> &imp_responce)
{
    assert(imp_responce.size()>0);
    if(kernelsize!=imp_responce.size())
    {
        return setKernel(imp_responce);
    }

    size_t nfft=cfg->nfft;
    size_t  n_imp_resp=imp_responce.size();

    assert((((int)nfft) - ((int)n_imp_resp) + 1)>=0);

    /*zero pad in the middle to left-rotate the impulse response
      This puts the scrap samples at the end of the inverse fft'd buffer */
    cfg->tmpbuf[0].r = imp_responce[ n_imp_resp - 1 ].real();
    cfg->tmpbuf[0].i = imp_responce[ n_imp_resp - 1 ].imag();
    for (size_t i=0;i<n_imp_resp - 1; ++i)
    {
        cfg->tmpbuf[ nfft - n_imp_resp + 1 + i ].r = imp_responce[ i ].real();
        cfg->tmpbuf[ nfft - n_imp_resp + 1 + i ].i = imp_responce[ i ].imag();
    }

    kiss_fft(cfg->fftcfg,cfg->tmpbuf,cfg->fir_freq_resp);

    /* TODO: this won't work for fixed point */
    float scale = 1.0 / ((float)cfg->nfft);

    for (size_t i=0; i < cfg->n_freq_bins; ++i ) {
    #ifdef USE_SIMD
        cfg->fir_freq_resp[i].r *= _mm_set1_ps(scale);
        cfg->fir_freq_resp[i].i *= _mm_set1_ps(scale);
    #else
        cfg->fir_freq_resp[i].r *= scale;
        cfg->fir_freq_resp[i].i *= scale;
    #endif
    }


    return nfft;


}

int QJFastFIRFilter::updateKernel(const QVector<kiss_fft_scalar> &imp_responce)
{
    assert(imp_responce.size()>0);
    if(kernelsize!=imp_responce.size())
    {
        return setKernel(imp_responce);
    }

    size_t nfft=cfg->nfft;
    size_t  n_imp_resp=imp_responce.size();

    assert((((int)nfft) - ((int)n_imp_resp) + 1)>=0);

    /*zero pad in the middle to left-rotate the impulse response
      This puts the scrap samples at the end of the inverse fft'd buffer */
    cfg->tmpbuf[0].r = imp_responce[ n_imp_resp - 1 ];
    cfg->tmpbuf[0].i = 0;
    for (size_t i=0;i<n_imp_resp - 1; ++i)
    {
        cfg->tmpbuf[ nfft - n_imp_resp + 1 + i ].r = imp_responce[ i ];
        cfg->tmpbuf[ nfft - n_imp_resp + 1 + i ].i = 0;
    }

    kiss_fft(cfg->fftcfg,cfg->tmpbuf,cfg->fir_freq_resp);

    /* TODO: this won't work for fixed point */
    float scale = 1.0 / ((float)cfg->nfft);

    for (size_t i=0; i < cfg->n_freq_bins; ++i ) {
    #ifdef USE_SIMD
        cfg->fir_freq_resp[i].r *= _mm_set1_ps(scale);
        cfg->fir_freq_resp[i].i *= _mm_set1_ps(scale);
    #else
        cfg->fir_freq_resp[i].r *= scale;
        cfg->fir_freq_resp[i].i *= scale;
    #endif
    }


    return nfft;


}

void QJFastFIRFilter::reset()
{
    kiss_fft_cpx asample;
    asample.r=0;
    asample.i=0;
    remainder.assign(nfft*2,asample);
    idx_inbuf=0;
    remainder_ptr=nfft;

    single_input_output_buf_ptr=0;
    single_input_output_buf.assign(nfft,asample);
}

void QJFastFIRFilter::Update(QVector<kiss_fft_cpx> &data)
{
    Update(data.data(), data.size());
}

void QJFastFIRFilter::Update(std::vector<kiss_fft_cpx> &data)
{
    Update(data.data(), data.size());
}

void QJFastFIRFilter::Update(kiss_fft_cpx *data,int Size)
{

    //ensure enough storage
    if((inbuf.size()-idx_inbuf)<(size_t)Size)
    {
        inbuf.resize(Size+nfft);
        outbuf.resize(Size+nfft);
    }

    //add data to storage
    memcpy ( inbuf.data()+idx_inbuf, data, sizeof(kiss_fft_cpx)*Size );
    size_t nread=Size;

    //fast fir of storage
    size_t nwrite=kiss_fastfir_complex(cfg, inbuf.data(), outbuf.data(),nread,&idx_inbuf);

    int currentwantednum=Size;
    int numfromremainder=std::min(currentwantednum,remainder_ptr);

    //return as much as posible from remainder buffer
    if(numfromremainder>0)
    {
        memcpy ( data, remainder.data(), sizeof(kiss_fft_cpx)*numfromremainder );

        currentwantednum-=numfromremainder;
        data+=numfromremainder;

        if(numfromremainder<remainder_ptr)
        {
            remainder_ptr-=numfromremainder;
            memcpy ( remainder.data(), remainder.data()+numfromremainder, sizeof(kiss_fft_cpx)*remainder_ptr );
        } else remainder_ptr=0;
    }

    //then return stuff from output buffer
    int numfromoutbuf=std::min(currentwantednum,(int)nwrite);
    if(numfromoutbuf>0)
    {
        memcpy ( data, outbuf.data(), sizeof(kiss_fft_cpx)*numfromoutbuf );
        currentwantednum-=numfromoutbuf;
        data+=numfromoutbuf;
    }

    //any left over is added to remainder buffer
    if(((size_t)numfromoutbuf<nwrite)&&(nwrite>0))
    {
        memcpy ( remainder.data()+remainder_ptr, outbuf.data()+numfromoutbuf, sizeof(kiss_fft_cpx)*(nwrite-numfromoutbuf) );
        remainder_ptr+=(nwrite-numfromoutbuf);
    }


    //if currentwantednum>0 then some items were not changed, this should not happen
    //we should anyways have enough to return but if we dont this happens. this should be avoided else a discontinuity of frames occurs. set remainder to zero and set remainder_ptr to nfft before running to avoid this
    if(currentwantednum>0)
    {
        qDebug()<<"Error: user wants "<<currentwantednum<<" more items from fir filter!";
        remainder_ptr+=currentwantednum;
    }

}

QJFastFIRFilter::~QJFastFIRFilter()
{
    free(cfg);
}

kiss_fft_cpx QJFastFIRFilter::Update_Single(kiss_fft_cpx signal)
{
    kiss_fft_cpx out_signal=single_input_output_buf[single_input_output_buf_ptr];
    single_input_output_buf[single_input_output_buf_ptr]=signal;
    single_input_output_buf_ptr++;single_input_output_buf_ptr%=single_input_output_buf.size();
    if(single_input_output_buf_ptr==0)Update(single_input_output_buf.data(),single_input_output_buf.size());
    return out_signal;
}

cpx_type QJFastFIRFilter::Update_Single(cpx_type signal)
{
    cpx_type out_signal=cpx_type(single_input_output_buf[single_input_output_buf_ptr].r,single_input_output_buf[single_input_output_buf_ptr].i);
    single_input_output_buf[single_input_output_buf_ptr].r=signal.real();
    single_input_output_buf[single_input_output_buf_ptr].i=signal.imag();
    single_input_output_buf_ptr++;single_input_output_buf_ptr%=single_input_output_buf.size();
    if(single_input_output_buf_ptr==0)Update(single_input_output_buf.data(),single_input_output_buf.size());
    return out_signal;
}

double QJFastFIRFilter::Update_Single(double signal)
{
    kiss_fft_cpx out_signal=single_input_output_buf[single_input_output_buf_ptr];
    kiss_fft_cpx tsig;tsig.i=0;tsig.r=signal;
    single_input_output_buf[single_input_output_buf_ptr]=tsig;
    single_input_output_buf_ptr++;single_input_output_buf_ptr%=single_input_output_buf.size();
    if(single_input_output_buf_ptr==0)Update(single_input_output_buf.data(),single_input_output_buf.size());
    return out_signal.r;
}

kiss_fft_cpx QJFastFIRFilter::Update_Single_c_out(double signal)
{


    kiss_fft_cpx out_signal=single_input_output_buf[single_input_output_buf_ptr];
    kiss_fft_cpx tsig;tsig.i=0;tsig.r=signal;
    single_input_output_buf[single_input_output_buf_ptr]=tsig;
    single_input_output_buf_ptr++;single_input_output_buf_ptr%=single_input_output_buf.size();
    if(single_input_output_buf_ptr==0)Update(single_input_output_buf.data(),single_input_output_buf.size());
    return out_signal;

}

cpx_type QJFastFIRFilter::Update_Single_c_out2(double signal)
{
    cpx_type out_signal=cpx_type(single_input_output_buf[single_input_output_buf_ptr].r,single_input_output_buf[single_input_output_buf_ptr].i);
    kiss_fft_cpx tsig;tsig.i=0;tsig.r=signal;
    single_input_output_buf[single_input_output_buf_ptr]=tsig;
    single_input_output_buf_ptr++;single_input_output_buf_ptr%=single_input_output_buf.size();
    if(single_input_output_buf_ptr==0)Update(single_input_output_buf.data(),single_input_output_buf.size());
    return out_signal;
}

//-----------

//real fast fir

QJFastFIRFilter_Real::QJFastFIRFilter_Real(QObject *parent) : QObject(parent)
{
    std::vector<kiss_fft_scalar> tvect;
    kernelsize=0;
    kiss_fft_scalar asample;
    asample=1;
    tvect.push_back(asample);
    nfft=2;
    cfg=kiss_fastfir_alloc_real(tvect.data(),tvect.size(),&nfft,0,0);
    reset();
}

QJFastFIRFilter_Real::~QJFastFIRFilter_Real()
{
    free(cfg);
}

void QJFastFIRFilter_Real::reset()
{
    kiss_fft_scalar asample;
    asample=0;
    remainder.assign(nfft*2,asample);
    idx_inbuf=0;
    remainder_ptr=nfft;


    single_input_output_buf_ptr=0;
    single_input_output_buf.assign(nfft,asample);

}

int QJFastFIRFilter_Real::setKernel(QVector<kiss_fft_scalar> imp_responce)
{
    int _nfft=imp_responce.size()*4;//rule of thumb
    _nfft=pow(2.0,(ceil(log2(_nfft))));
    return setKernel(imp_responce,_nfft);
}

int QJFastFIRFilter_Real::setKernel(QVector<kiss_fft_scalar> imp_responce,int _nfft)
{
    if(!imp_responce.size())return nfft;
    free(cfg);
    _nfft=pow(2.0,(ceil(log2(_nfft))));
    nfft=_nfft;
    cfg=kiss_fastfir_alloc_real(imp_responce.data(),imp_responce.size(),&nfft,0,0);
    reset();
    kernelsize=imp_responce.size();
    return nfft;
}

void QJFastFIRFilter_Real::Update(std::vector<kiss_fft_scalar> &data)
{
    Update(data.data(), data.size());
}

void QJFastFIRFilter_Real::Update(kiss_fft_scalar *data,int Size)
{

    //ensure enough storage
    if((inbuf.size()-idx_inbuf)<(size_t)Size)
    {
        inbuf.resize(Size+nfft);
        outbuf.resize(Size+nfft);
    }

    //add data to storage
    memcpy ( inbuf.data()+idx_inbuf, data, sizeof(kiss_fft_scalar)*Size );
    size_t nread=Size;

    //fast fir of storage
    size_t nwrite=kiss_fastfir_real(cfg, inbuf.data(), outbuf.data(),nread,&idx_inbuf);

    int currentwantednum=Size;
    int numfromremainder=std::min(currentwantednum,remainder_ptr);

    //return as much as posible from remainder buffer
    if(numfromremainder>0)
    {
        memcpy ( data, remainder.data(), sizeof(kiss_fft_scalar)*numfromremainder );

        currentwantednum-=numfromremainder;
        data+=numfromremainder;

        if(numfromremainder<remainder_ptr)
        {
            remainder_ptr-=numfromremainder;
            memcpy ( remainder.data(), remainder.data()+numfromremainder, sizeof(kiss_fft_scalar)*remainder_ptr );
        } else remainder_ptr=0;
    }

    //then return stuff from output buffer
    int numfromoutbuf=std::min(currentwantednum,(int)nwrite);
    if(numfromoutbuf>0)
    {
        memcpy ( data, outbuf.data(), sizeof(kiss_fft_scalar)*numfromoutbuf );
        currentwantednum-=numfromoutbuf;
        data+=numfromoutbuf;
    }

    //any left over is added to remainder buffer
    if(((size_t)numfromoutbuf<nwrite)&&(nwrite>0))
    {
        memcpy ( remainder.data()+remainder_ptr, outbuf.data()+numfromoutbuf, sizeof(kiss_fft_scalar)*(nwrite-numfromoutbuf) );
        remainder_ptr+=(nwrite-numfromoutbuf);
    }


    //if currentwantednum>0 then some items were not changed, this should not happen
    //we should anyways have enough to return but if we dont this happens. this should be avoided else a discontinuity of frames occurs. set remainder to zero and set remainder_ptr to nfft before running to avoid this
    if(currentwantednum>0)
    {
        qDebug()<<"Error: user wants "<<currentwantednum<<" more items from fir filter!";
        remainder_ptr+=currentwantednum;
    }

}

double QJFastFIRFilter_Real::Update_Single(kiss_fft_scalar signal)
{
    kiss_fft_scalar out_signal=single_input_output_buf[single_input_output_buf_ptr];
    single_input_output_buf[single_input_output_buf_ptr]=signal;
    single_input_output_buf_ptr++;if(single_input_output_buf_ptr>=single_input_output_buf.size())single_input_output_buf_ptr=0;
    if(single_input_output_buf_ptr==0)Update(single_input_output_buf.data(),single_input_output_buf.size());
    return out_signal;
}

int QJFastFIRFilter_Real::updateKernel(const QVector<kiss_fft_scalar> &imp_responce)
{
    assert(imp_responce.size()>0);
    if(kernelsize!=imp_responce.size())
    {
        return setKernel(imp_responce);
    }

    size_t nfft=cfg->nfft;
    size_t  n_imp_resp=imp_responce.size();

    assert((((int)nfft) - ((int)n_imp_resp) + 1)>=0);

    /*zero pad in the middle to left-rotate the impulse response
      This puts the scrap samples at the end of the inverse fft'd buffer */
    cfg->tmpbuf[0] = imp_responce[ n_imp_resp - 1 ];
    for (size_t i=0;i<n_imp_resp - 1; ++i)
    {
        cfg->tmpbuf[ nfft - n_imp_resp + 1 + i ] = imp_responce[ i ];
    }

    kiss_fftr(cfg->fftcfg,cfg->tmpbuf,cfg->fir_freq_resp);

    /* TODO: this won't work for fixed point */
    float scale = 1.0 / ((float)cfg->nfft);

    for (size_t i=0; i < cfg->n_freq_bins; ++i ) {
    #ifdef USE_SIMD
        cfg->fir_freq_resp[i].r *= _mm_set1_ps(scale);
        cfg->fir_freq_resp[i].i *= _mm_set1_ps(scale);
    #else
        cfg->fir_freq_resp[i].r *= scale;
        cfg->fir_freq_resp[i].i *= scale;
    #endif
    }


    return nfft;


}

//------

//---fast Hilbert filter (actually a analytical function creator so it is a hilbert transform multiplied by imaginary one and this is added to the real input)
QJHilbertFilter::QJHilbertFilter(QObject *parent)  : QJFastFIRFilter(parent)
{
    setSize(2048);
}

void QJHilbertFilter::setSize(int N)
{
    N=pow(2.0,(ceil(log2(N))));
    assert(N>1);

    kernel.clear();

    cpx_type asample;
    for(int i=0;i<N;i++)
    {

        if(i==N/2)
        {
            //asample.i=0;
            //asample.r=-1;
            asample=cpx_type(-1,0);
            kernel.push_back(asample);
            continue;
        }

        if((i%2)==0)
        {
            //asample.i=0;
            //asample.r=0;
            asample=cpx_type(0,0);
            kernel.push_back(asample);
            continue;
        }

        //asample.r=0;
        //asample.i=(2.0/((double)N))/(std::tan(M_PI*(((double)i)/((double)N)-0.5)));
        asample=cpx_type(0,(2.0/((double)N))/(std::tan(M_PI*(((double)i)/((double)N)-0.5))));
        kernel.push_back(asample);

    }

    //flip real imaginary. for some reason I got the kernel calulations back to front
    // i think that was due to not being able to spin the phasors the other way
    for(int k=0;k<kernel.size();k++)kernel[k]=cpx_type(kernel[k].imag(),kernel[k].real());

    setKernel(kernel);
}

QVector<cpx_type> QJHilbertFilter::getKernel()
{
    return kernel;
}

void QJHilbertFilter::setSecondFilterKernel(QVector<kiss_fft_scalar> imp_responce)
{
    //convolute the two kernels
    QVector<cpx_type> imp_responce_combined;
    imp_responce_combined.resize(imp_responce.size()+kernel.size()-1);
    QVector<cpx_type> buff;
    buff.fill(0,kernel.size());
    for(int i=0;i<imp_responce_combined.size();i++)
    {
        for(int k=buff.size()-2;k>=0;k--)buff[k+1]=buff[k];
        if(i<imp_responce.size())buff[0]=imp_responce[i];
         else buff[0]=0;
        cpx_type val=0;
        for(int j=0;j<kernel.size();j++)
        {
            val+=buff[j]*kernel[j];
        }
        imp_responce_combined[i]=val;
    }
    updateKernel(imp_responce_combined);
}

void QJHilbertFilter::setSecondFilterKernel(QVector<cpx_type> imp_responce)
{
    //convolute the two kernels
    QVector<cpx_type> imp_responce_combined;
    imp_responce_combined.resize(imp_responce.size()+kernel.size()-1);
    QVector<cpx_type> buff;
    buff.fill(0,kernel.size());
    for(int i=0;i<imp_responce_combined.size();i++)
    {
        for(int k=buff.size()-2;k>=0;k--)buff[k+1]=buff[k];
        if(i<imp_responce.size())buff[0]=imp_responce[i];
         else buff[0]=0;
        cpx_type val=0;
        for(int j=0;j<kernel.size();j++)
        {
            val+=buff[j]*kernel[j];
        }

        imp_responce_combined[i]=val;
    }
    updateKernel(imp_responce_combined);
}


//---

