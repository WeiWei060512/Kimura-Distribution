#!/usr/bin/env python
import sys
import csv
import pandas as pd
#import glob
import numpy as np
import random
#import timeit

##Calculating numerical values of hypergeometric function, F(i-1, i+2, 2, x)
##The algorithm for calculating this function is based on Hoang-Binh 2005.

def generate_hypergeometric_function(id, x):
    a = 1-id
    b = id +2
    c = 2

    f0 = 1.0
    f1 = 1-(b*x/c)

    if(id == 1):
        return f0

    if(id == 2):
        return f1
    
    f = 0.0
    if(id > 2):
        for j in range(2, id):
            aa = 1-j
            f = (aa*(1-x)*(f1-f0)+(aa+b*x-c)*f1)/(aa-c)
            f0 = f1
            f1 = f
        
        return f    
        
    
#Calculate phi function for each x value
#This function was called by function finding_converge_id.
#The calculation of infinite phi function will be stopped when the numerical value of phi(x) at the ith term is less than 0.0001. 
#[At the ith term, phi(x) < 0.0001 is equivalent to (summation of phi(x) at the (i+1)th - summation of phi(x) at the ith term) < 0.0001].   

def test_through_generate_phi_function(id, p, x, b_ne):
    q = 1-p
    sum_phi = 0.0       
    phi=list()

    for i in range(1, id+1):
        
        #1st component: pqi(i+1)(2i+1)
        pq_coeff = p*q*i*(i+1)*((2*i)+1)
        
        #2nd component: F(1-i,i+2,2,p)F(1-i,i+2,2,x)    
        pp_hypgeo = generate_hypergeometric_function(i, p)
        xx_hypgeo = generate_hypergeometric_function(i, x)       
        
        #3rd component in the form of b parameter: b^[i(i+1)/2]
        b_coeff = int(i*(i+1)/2)
        b_exp = 1.0
        for j in range(1, b_coeff+1):
            b_exp *= b_ne

        #phi(x,b): put all components together
        phi_list = pq_coeff*pp_hypgeo*xx_hypgeo*b_exp
        phi.append(phi_list)
        
        #test criterium for truncating the infinite series of Kimura's probability density function
        if((abs(phi_list) < 0.00001) and (i > 1)):
            use_id = i
            return use_id
            break
        sum_phi += phi_list

        
    use_id = id
    
    return use_id


#Getting the maximum number of the iterations for calculating Kimura distribution
#This number is obtained from comparing among 10 ith terms provided by calculating phi function based on 10 x values. 
#Each x-vlaue is varied from 0 to 1.0 which an increment is equal to 0.1. 

def finding_converge_id(terms, p, b_ne):
    converge_id=list() 
    for k in range(0, 110, 10):
        x = k/100.00
        h = k/10
        converge_id_list=test_through_generate_phi_function(terms, p, x, b_ne)
        converge_id.append(converge_id_list) 
    sorted(converge_id)
    use_id = converge_id [9]
        
    return use_id;


#Calculating probability of x = 0 (loss) or x = 1 (fixed) based on Kimura distribution [ref 23] 

def generate_prob_loss_and_fixed(converge_id, pq, b_ne):
    sum_fixed = 0.0
    for i in range(1, converge_id+1):
    
        #1st component: pq(2i+1)(-1)^[i]
        if (i %2 == 0):
            lf_equation_coeff = 1*((2*i)+1)
        else:
            lf_equation_coeff = -1*((2*i)+1) 
        lf_pq_coeff = pq*(1-pq)*lf_equation_coeff
      
        #2nd component: F(1-i,i+2,2,p)  
        lf_pq_hypgeo = generate_hypergeometric_function(i, pq)
                
        #3rd component: b^[i(i+1)/2]
        b_coeff = i*(i+1)/2
        b_coeff = int(b_coeff)
        b_exp = 1.0 
        for j in range(1, b_coeff+1):
            b_exp *= b_ne
        
        #f(1,b): put all components together    
        fixed = lf_pq_coeff*lf_pq_hypgeo*b_exp
        sum_fixed += fixed
    
    fixed_prob = pq + sum_fixed 
    return fixed_prob   
    
#Calculating phi function based on Kimura distribution [ref 23]
def generate_phi_function(converge_id, p, x, b_ne):
    q = 1-p
    sum_phi = 0.0
    for i in range(1, converge_id+1): 
        #1st component: pqi(i+1)(2i+1)
        pq_coeff = p*q*i*(i+1)*((2*i)+1)
        #2nd component: F(1-i,i+2,2,p)F(1-i,i+2,2,x)    
        pp_hypgeo = generate_hypergeometric_function(i, p)
        xx_hypgeo = generate_hypergeometric_function(i, x)
        #3rd component in the form of b parameter: b^[i(i+1)/2]
        b_coeff = i*(i+1)/2
        b_coeff = int(b_coeff)
        b_exp = 1.0
        for j in range(1, b_coeff+1):
            b_exp *= b_ne
        #phi(x,b): put all components together  
        phi = pq_coeff*pp_hypgeo*xx_hypgeo*b_exp
        sum_phi += phi
        
    return sum_phi

def read_input_data(inputfile, outputfile, sample_size, limit_size):
    het_score = list()
    het_level = list()
    with open(inputfile) as csvfile:
        reader = csv.reader(csvfile, delimiter='\t')
        for row in reader:
            het_score.append(float(row[0]))
    for u in range(0, limit_size + 1):
        het_level.append(0)
    for cscore_int in range(0, limit_size + 1):
        cscore = cscore_int / (limit_size * 1.00)
        for pscore in range(0, sample_size):
            if ((het_score[pscore] <= cscore) and (het_score[pscore] >= 0)):
                het_level[cscore_int] += 1;
    with open(outputfile, 'w') as csvfile:
        csv_writer = csv.writer(csvfile, delimiter='\t')
        for w in range (0, limit_size + 1):
            csv_writer.writerow([het_level[w] / (sample_size * 1.00)])
    csvfile.close

def read_input_data_2(sim_data, limit_max, sample_size):
    het_level = list()
    het_level_sim = list()
    for u in range(0, limit_max + 1):
        het_level.append(0)
    for chet_int in range(0, limit_max + 1):
        chet = chet_int / (limit_max * 1.00)
        for sd in range(0, sample_size):
            if ((sim_data[sd] <= chet) and (sim_data[sd] >= 0)):
                het_level[chet_int] += 1
    for v in range(0, limit_max + 1):
        het_level_sim.append(het_level[v] / (sample_size * 1.00))
    return het_level_sim

#Applying trapzoid rule to integrate phi function
def trapezoid(converge_id, aa, bb, n, p, b_ne):
    del_del = (bb-aa)/n

    func_a = generate_phi_function(converge_id, p, aa, b_ne)
    func_b = generate_phi_function(converge_id, p, bb, b_ne)
        
    sum1 = (0.5*del_del)*(func_a + func_b)
    sum2 = 0.00

    for j in range(1, n):
        y = aa + (j*del_del)
        func_y = generate_phi_function(converge_id,p,y,b_ne)
        sum2 += func_y
        
    s = sum1 + (del_del*sum2)
    return s

#Calculating numerical values of cumulative probability distribution function (cdf) based on Kimura distribution [ref 23]
#The result will be written into the file. The name of this file was detrmined in char theoryfile.
#Specifically, the purpose of this function is generating cdf of continuous variables caclulated based on Kimura probability density function which will be applied to generate simulated data sets. 

def calculate_continuous_cdf(output, converge_id, p, b_ne, sub_integ, max_data):
    prob_cdf=list() 
    prob_integrated=list()
    prob_integrated.append(0)
    q = 1 - p
    prob_cdf.append(generate_prob_loss_and_fixed(converge_id, q, b_ne))
    f = open(output,'w')
    f.write("bin\tP_kimura\n") 
    f.write("{0:10d}\t{1:10.6f}\n".format(0, prob_cdf[0]))
    
    dividen = max_data/1.00
    aa_int = 0
    aa = aa_int/dividen

    for bb_int in range(1, max_data):
        bb = bb_int/dividen
        prob_integrated.append(trapezoid(converge_id, aa, bb, sub_integ, p, b_ne))
        prob_cdf.append(prob_cdf[0] + prob_integrated[bb_int])
        f.write("{0:10d}\t{1:10.6f}\n".format(bb_int, prob_cdf[bb_int]))

    prob_hetlast_integrated = trapezoid(converge_id, 0.0, 1.0, sub_integ, p, b_ne)
    prob_hetlast = prob_cdf[0] + prob_hetlast_integrated
    
    prob_fixed = generate_prob_loss_and_fixed(converge_id, p, b_ne)
        
    prob_cdf.append(prob_hetlast + prob_fixed)
    f.write("{0:10d}\t{1:10.6f}\n".format(max_data, prob_cdf[max_data]))
    with open(output + '2', 'w') as csvfile:
        csv_writer = csv.writer(csvfile, delimiter='\t')
        csv_writer.writerows(map(lambda x: [x], prob_cdf))
    csvfile.close()
    return prob_cdf
    f.close(out_file)

#Calculating numerical values of cumulative probability distribution function (cdf) based on Kimura distribution [ref 23]
#The result will be written into the file. The name of this file was detrmined in char theoryfile.
#Specifically, the purpose of this function is generating the fit cumulative Kimura probability distribution for drawing a histogram
#Specifically, the purpose of this function is generating cdf of discrete variables caclulated based on Kimura probability density function which will be applied to compare to simulated data sets and observed data set.   

def calculate_discrete_cdf(output, converge_id, p, b_ne, sub_integ, no_binned):
    prob_cdf=list() 
    prob_integrated=list()
    q = 1-p
    prob_cdf.append(generate_prob_loss_and_fixed(converge_id, q, b_ne))
    f = open(output,'w')
    f.write("bin\tP_kimura\n") 
    f.write("{0:5.4f}\t{1:10.6f}\n".format(0.00, prob_cdf[0]))
    aa_int = 0
    aa = aa_int/100.00
    for bb_int in range(10, 100, 10):
        bb = bb_int/100.00
        h = int((bb_int/10)-1)
        prob_integrated.append(trapezoid(converge_id, aa, bb, sub_integ, p, b_ne))
        prob_cdf.append(prob_cdf[0] + prob_integrated[h])
        f.write("{0:5.4f}\t{1:10.6f}\n".format(bb, prob_cdf[h+1]))

    prob_all_het = trapezoid(converge_id, 0.00, 1.00, sub_integ, p, b_ne) 
    prob_fixed = generate_prob_loss_and_fixed(converge_id, p, b_ne)     
    prob_cdf [9] = prob_cdf[0] + prob_all_het + prob_fixed
    f.write("{0:5.4f}\t{1:10.6f}\n".format(1.0, prob_cdf[9]))
      
    return prob_cdf
    f.close() 


#Calculating maximum values of the difference between numerical values of observed cdf and kimura cdf 

def Calculating_max_distance1(observed_comp, kimura_comp, output):
    df_merged_comp = pd.DataFrame({'observed_dis':observed_comp, 'kimura_dis': kimura_comp })
    df_merged_comp['D'] = df_merged_comp['kimura_dis'] - df_merged_comp['observed_dis']
    df_merged_comp['abs_D'] = abs(df_merged_comp['kimura_dis'] - df_merged_comp['observed_dis'])
    Max_D = max(df_merged_comp['abs_D'])
    df_merged_comp.to_csv(output, sep='\t', index = False, na_rep='')
    return Max_D

#Calculating maximum values of the difference between numerical values of simulated cdf and kimura cdf 

def Calculating_max_distance2(simu_comp, kimura_comp):
    df_merged_comp = pd.DataFrame({'simu_dis': simu_comp, 'kimura_dis': kimura_comp})
    df_merged_comp['D'] = df_merged_comp['kimura_dis'] - df_merged_comp['simu_dis']
    df_merged_comp['abs_D'] = abs(df_merged_comp['kimura_dis'] - df_merged_comp['simu_dis'])
    Max_D = max(df_merged_comp['abs_D'])
    return Max_D

def sim_ks_test(sim_cdf, comp_phi_cdf, limit_size):
    diff1 = list()
    for sd in range(0, limit_size + 1):
        diff1.append(0.00)
    max_diff1 = 0.00

    for sc in range(0, limit_size + 1):
        #calculate the difference between numerical values of observed cdf and expected cdf
        #diff1.append(abs(comp_phi_cdf[sc] - sim_cdf[sc]))
        diff1[sc] = abs(comp_phi_cdf[sc] - sim_cdf[sc])
        if (diff1[sc] >= max_diff1):
            max_diff1 = diff1[sc]    
    return max_diff1

def linear_search(array, key, size):
    if(key < array[0]):
        return 0
    else:
        for n in range(1, size):
            m = n + 1
            lin_est = abs(key-array[n])

            if (abs(array[m] - array[n]) >= lin_est):
                x_value = n + lin_est
                return x_value
    return -1

def calculate_simulated_cdf(sim_data, sample_size, run, no_binned):
    het_level = list()
    het_level_sim = list()
    for u in range(0, limit_max+1):
        het_level.append(0)

    for chet_int in range(0, limit_max+1):
        chet = chet_int/((limit_max)*1.00)
        for sd in range(0, sample_size):
             if((sim_data[sd] <= chet) and (sim_data[sd] >= 0)):
                het_level[chet_int] += 1
    for v in range(0, limit_max+1):
        het_level_sim.append(het_level[v-1]/(sample_size*1.00))

    return het_level_sim

def main():
    if len(sys.argv) < 2:
        print("You must set input raw data without header please!!!\n")
    elif len(sys.argv) < 3:
        print("You must set output prefix with full path!!!\n")
    elif len(sys.argv) == 3:
        print("Well done, input and output files are tracked, check your input including one column (float) & without header, outputs with full path!!!") 
    
    input_raw = sys.argv[1]
    output_pre = sys.argv[2]
    output_summary = output_pre + '.summary.txt' 
    output_binned1 = output_pre + '.observed.barchart.csv'
    output_binned2 = output_pre + '.kimura.simu.csv'   ## kimura cumuplot input
    output_binned3 = output_pre + '.kimura.comp.csv'
    output_binned4 = output_pre + '.kimura.barchart.csv'
    output_binned5 = output_pre + '.observed.comp.csv'
    output_binned6 = output_pre + '.simu.cdf.txt'
    output_binned7 = output_pre + '.observed.cumuplot.csv'

    output_monte_carlo1000 = output_pre + '.monte_carlo1000.txt'
    output_barchart = output_pre + '.barchart.pdf'
    output_cumuplot = output_pre + '.cumuDist.pdf'

    #calculate input parameters: var = heteroplasmy variance, p = p initial heteroplasmy, b_ne = b parameter

    df=pd.read_csv(input_raw, names=['HF'], sep="\t")
    df.sort_values(by=['HF'], inplace=True)
    p=df['HF'].mean()
    var=np.var(df['HF'])
    sample_size = len(df)
    b_ne = ((p*(1-p)) - var)/ (p*(1-p))

#   var = 0.0143;
#    p = 0.1264;
#    sample_size = 82
#    b_ne = ((p*(1-p)) - var)/ (p*(1-p))

    #Calculate the maxinum number of iterations needed to be used in calculating Kimura distribution
    terms = 10;
    converge_id = finding_converge_id(terms, p, b_ne)
    
    #print and write to output1
    sys.stdout=open(output_summary,"w")
    print("Parameters: p is {0:5.4f} with variance = {1:5.4f}, b is {2:5.4f}\n".format(p, var, b_ne))
    print ("Converge_id [@sum_phi < 0.0001]: {0:d}\n".format(converge_id))

    #Resolution of Kimura distribution and integrated Kimura distribution
    sub_integ = 200 #The number of subdivisions needed for integrating Kimura distribution 
    max_data = 10 #This number determines resolution of calculating Kimura distribution for generating simulated data 
    no_binned = 10 #Number of binned data in a data set
    limit_size = max_data * 10 #This number determines resolution of calculating Kimura distribution for                                                            comparing to the observed data and calculating cdf of simulatied data 

    #Calculate Kimura distribution 

    #1. cdf of continuous variables: for generating simulated data set
    try:
        phi_cdf = list()
        with open(output_binned2 + '2') as csvfile:
            reader = csv.reader(csvfile, delimiter='\t')
            for row in reader:
                phi_cdf.append(float(row[0]))
        csvfile.close()
    except IOError:
        phi_cdf = calculate_continuous_cdf(output_binned2, converge_id, p, b_ne, sub_integ, max_data)

    #2.  cdf of continuous variables: for compraring data set 
    
    try:
        comp_phi_cdf = list()
        with open(output_binned3 + '2') as csvfile:
            reader = csv.reader(csvfile, delimiter='\t')
            for row in reader:
                comp_phi_cdf.append(float(row[0]))
        csvfile.close()
    except IOError:
        comp_phi_cdf = calculate_continuous_cdf(output_binned3, converge_id, p, b_ne, sub_integ, limit_size)
 
    sample_cdf = list()
    read_input_data(input_raw, output_binned5, sample_size, limit_size);
    with open(output_binned5) as csvfile:
        reader = csv.reader(csvfile, delimiter='\t')
        for row in reader:
            sample_cdf.append(float(row[0]))
    sample_D = sim_ks_test(sample_cdf, comp_phi_cdf, limit_size);
    print("sample_D: ", sample_D)

    #Generate simulated dataset
    max_D = list()
    random.seed(-942320)

    f_monto = open(output_monte_carlo1000,'w')
    f_monto.write("Sample:")
    f_monto.write("{0:0.6f}\n".format(sample_D)) 

    #The number of simulated dataset in each run is determined by max_run.              
    max_run = 1000
    higher_t_sample = 0
    low_or_eq = 0
        
    #Generate simulated dataset of each repeated run with the number of dataset is equal to max_run
    with open(output_pre + ".monte_carlo.txt", 'w') as csvfile:
        csv_writer = csv.writer(csvfile, delimiter='\t')
        for run in range(0, max_run):
            art_mea = list()
            art_data = list()
            for numr in range(0, sample_size):
                rand_prob = random.random()
                art_data.append(linear_search(phi_cdf, rand_prob, max_data))
                art_mea.append(art_data[numr] / (max_data * 1.0))
            sim_cdf = read_input_data_2(art_mea, limit_size, sample_size)
            max_D.append(sim_ks_test(sim_cdf, comp_phi_cdf, limit_size))
            csv_writer.writerow([run, max_D[run], sample_D, max_D[run] - sample_D])

            #comparing maximum differences obtained from simulated data to the observed data
            if((max_D[run] - sample_D) >= 0.0000):
                f_monto.write("{0:d}\t{1:7.6f}\t{2:7.6f}\t{3:7.6f}\n".format(run, max_D[run], sample_D, max_D[run]-sample_D))
                higher_t_sample += 1    
            else:
                f_monto.write("{0:d}\t{1:7.6f}\t{2:7.6f}\t{3:7.6f}\n".format(run, max_D[run], sample_D, max_D[run]-sample_D))
                low_or_eq += 1
        csv_writer.writerow(["SUM", higher_t_sample, low_or_eq, higher_t_sample + low_or_eq])
    f_monto.write("Monte Carlo:\t{0:d}\t{1:d}\t{2:d}\n".format(higher_t_sample, low_or_eq, higher_t_sample+low_or_eq))


    #p-value of the KS test is eqaul the proportion of the maximum differences obtained from simulated data that are higher than the one obtained from observed data 
    f_monto.write("P-value: ")
    f_monto.write("{0:0.6f}\n".format(higher_t_sample/(max_run * 1.00)))
    print("P-value: ", higher_t_sample/(max_run * 1.00))
    f_monto.close()
   
 

if __name__ == '__main__':
    main() 

