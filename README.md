# NCTU Advanced Communication System Simulations and Experiments
## Lab 1 MATLAB introduction



## [Lab 2 Signals and the Fourier transform](https://github.com/EthanLiao/NCTU-Advanced-Communication-System-Simulations-and-Experiments/tree/master/Lab_2)

- Practice 1 to 5
    - Practice 1-2 : sinusoidal wave with different frequency
    - Practice 3 : Padding with Signal
    - Practice 4 : Rectangular Wave
    - Practice 5 : Triangular Wave With AWGN

    ![](https://i.imgur.com/iBafg1k.png)

    - HW1 : Sinusoidal Wave with High Sampling Frequency

    ![](https://i.imgur.com/RkoYbHM.png)

    - HW2 : Use Window to filter triangular wave noise

    ![](https://i.imgur.com/JXhktbl.png)

## [Lab 3 LTI systems and their responses](https://github.com/EthanLiao/NCTU-Advanced-Communication-System-Simulations-and-Experiments/tree/master/Lab_3)

- practice 1 and homework : designed IIR and FIR filter
    ![](https://i.imgur.com/rLZJUSD.png)
    ![](https://i.imgur.com/tqvNTZs.png)

## [Lab 4 Analog modulation ](https://github.com/EthanLiao/NCTU-Advanced-Communication-System-Simulations-and-Experiments/tree/master/Lab_4)

- practice 1 :
    ![](https://i.imgur.com/P8NwRNy.png)

    ![](https://i.imgur.com/XXeWhKT.png)

    #### practice 1 fail : carrier frequency = 1/64
    ![](https://i.imgur.com/Q97J7uj.png)

    #### practice 1 success : carrier frequency = 1/4
    ![](https://i.imgur.com/8F8ombT.png)

- practice 2 : 
    ![](https://i.imgur.com/CtFhaAw.png)
    ![](https://i.imgur.com/qjRKUOL.png)

- homework : integrate practice one and two then considering delay in channel
    ![](https://i.imgur.com/m6eQ0KR.png)


## [Lab 5 Digital modulation ](https://github.com/EthanLiao/NCTU-Advanced-Communication-System-Simulations-and-Experiments/tree/master/Lab_5)
- practice 1 : Generate 8-PAM and 16-QAM mapper
   
    - 8 -PAM mapper :
    ![](https://i.imgur.com/ZvrdBMN.png)

    ![](https://i.imgur.com/POWisGk.png)

    - 16-QAM mapper :
    ![](https://i.imgur.com/CcCoPvG.png)
    ![](https://i.imgur.com/g00baJs.png)

    ![](https://i.imgur.com/bhFoNtm.png)

- practice 2 : Generate a 16-QAM sequence with SNR 10dB
    ![](https://i.imgur.com/kYUB4do.png)

- practice 3 : Conduct a detection and calculate symbol error rate
    ![](https://i.imgur.com/BbOsNgO.png)

- HW : Symulate the symbol error rate of the 16-QAM with 5dB,10dB,15dB and then calculate the theoratical SER
    ![](https://i.imgur.com/kXertEU.png)

## [Lab 6 Sampling and rate conversion ](https://github.com/EthanLiao/NCTU-Advanced-Communication-System-Simulations-and-Experiments/tree/master/Lab_6)

- practice 1 : Generate a down sampling and up sampling signal and observe its spectrum
    ![](https://i.imgur.com/odXx2SN.png)

- practice 2 : up sample a signal and use IIR and FIR filter to recover it
    ![](https://i.imgur.com/FTcu9Pn.png)

- practice 3 : up then down sample a signal and use IIR and FIR filter to recover it
    - cut-off = 0.4 
    ![](https://i.imgur.com/KLYXMb2.png)
    ![](https://i.imgur.com/xJPz7H0.png)


- HW : dwon then up sample a signal and use IIR and FIR filter to recover it
    - cut-off = 0.5
    ![](https://i.imgur.com/A2T9ZTx.png)
    
    ![](https://i.imgur.com/lhAbvGg.png)

    ![](https://i.imgur.com/2dbn7BV.png)

    ![](https://i.imgur.com/Ft0X6DN.png)
    
    ![](https://i.imgur.com/5bXIEDG.png)

## [Lab 7 Transmit filtering/up conversion I ](https://github.com/EthanLiao/NCTU-Advanced-Communication-System-Simulations-and-Experiments/tree/master/Lab_7)
   - practice 1 : plot raised cosine and square root cosine
   ![](https://i.imgur.com/0mqpc5a.png)

   - practice 2 : Conduct srrc pulse shapping
    ![](https://i.imgur.com/QhMFxOl.png)
    ![](https://i.imgur.com/lvBzA6d.png)
    
        ![](https://i.imgur.com/9GN31pO.png)

   - practice 3 : Conduct filter pulse shapping
       - cut-off frequency = 0.4
       ![](https://i.imgur.com/m2DLyX4.png)
       ![](https://i.imgur.com/4OQNi9S.png)

   - HW : QPSK-UP conversion:
       ![](https://i.imgur.com/aApqECj.png)
        ![](https://i.imgur.com/XJq1tgf.png)

    
## [Lab 8 Transmit filtering/up conversion II](https://github.com/EthanLiao/NCTU-Advanced-Communication-System-Simulations-and-Experiments/tree/master/Lab_8)
    
- practice 1 : use SRRC as DAC and DMA filter


![](https://i.imgur.com/evJlbEg.png)
![](https://i.imgur.com/dR4AzAA.png)

- practice 2 :
![](https://i.imgur.com/z2EylfH.png)
![](https://i.imgur.com/FQmiDdP.png)


- practice 3 : design a IIR DMA filter then replace the previous problem DMA filter
    - cut-off : 
        - pass band :0.125 ; stop band : 0.4375
    ![](https://i.imgur.com/uK8HQQt.png)

- HW : 
    - cut-off : 
        - pass band :0.6 ; stop band : 0.825 ; Apass : 1 ; Astop: 60
![](https://i.imgur.com/Z5AAMFP.png)

## [Lab 9 Receive filtering/down conversion](https://github.com/EthanLiao/NCTU-Advanced-Communication-System-Simulations-and-Experiments/tree/master/Lab_9)


## [Lab 10 RF impairments](https://github.com/EthanLiao/NCTU-Advanced-Communication-System-Simulations-and-Experiments/tree/master/Lab_10)


## [Lab 11 Equalization ](https://github.com/EthanLiao/NCTU-Advanced-Communication-System-Simulations-and-Experiments/tree/master/Lab_11)


## [Lab 12 Constant envelop modulation ](https://github.com/EthanLiao/NCTU-Advanced-Communication-System-Simulations-and-Experiments/tree/master/Lab_12)


## [Lab 13 MIMO transmission](https://github.com/EthanLiao/NCTU-Advanced-Communication-System-Simulations-and-Experiments/tree/master/Lab_13)


## [Lab 14 Fixed-point implementation]()
## [Lab 15 Testing]()
