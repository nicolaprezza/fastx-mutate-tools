package fastx_mutate_tools;

import java.io.*;

public class FastxMutateTools {

	static public enum program_mode{
		snp_sim,indel_sim,bs_sim
	}
	
	static String input="";
	static String output="";
	static String output_meth="";
	static double snp=0;
	static double bs_fwd=0;
	static double bs_rev=0;
	static double cr=0.8;
	static int start_coord=0;
	static double open=0.005;
	static double extend=0.2;
	
	static program_mode mode;
	
	static BufferedReader in_fastx=null;//input
	static BufferedWriter out_fastx=null;
	static BufferedWriter out_meth_annotation=null;
	
	static void initBuffers(){
		
		if(input.compareTo("")==0){
						
			in_fastx = new BufferedReader(new InputStreamReader(System.in));
			
		}else{
		
			try{
				in_fastx = new BufferedReader(new InputStreamReader(new FileInputStream(input)));
			}catch(FileNotFoundException e){
				out("Error: file " + input + " not found. Usage:");
				help();
			}
			
		}
		
		if(output.compareTo("")==0){
			
			out_fastx = new BufferedWriter(new OutputStreamWriter(System.out));

		}else{
			
			try{
				out_fastx = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(output)));
			}catch(FileNotFoundException e){
				out("Error: file " + output + " not found. Usage:");
				help();
			}
			
		}
		
		if(output_meth.compareTo("")!=0){
			
			try{
				out_meth_annotation = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(output_meth)));
			}catch(FileNotFoundException e){
				out("Error: file " + output_meth + " not found. Usage:");
				help();
			}
			
		}
		
	}
	
	static void out(String s){
		System.out.println(s);
	}
	
	static int parse(String[] args, int i){
		
		if(args[i].compareTo("--open")==0){
			i++;
			open = Double.parseDouble(args[i]);
			i++;			
		}else if (args[i].compareTo("--extend")==0){
			i++;
			extend = Double.parseDouble(args[i]);
			i++;
		}else if(args[i].compareTo("--snp")==0){
			i++;
			snp = Double.parseDouble(args[i]);
			i++;			
		}else if (args[i].compareTo("--1-based")==0){
			start_coord=1;
			i++;
		}else if (args[i].compareTo("--bs-fwd")==0){
			i++;
			bs_fwd = Double.parseDouble(args[i]);
			i++;
		}else if (args[i].compareTo("--bs-rev")==0){
			i++;
			bs_rev = Double.parseDouble(args[i]);
			i++;			
		}else if (args[i].compareTo("--output")==0){
			i++;
			output = args[i];
			i++;
		}else if (args[i].compareTo("--input")==0){
			i++;
			input = args[i];
			i++;
		}else if (args[i].compareTo("--cr")==0){
			i++;
			cr = Double.parseDouble(args[i]);
			i++;
		}else if (args[i].compareTo("--output-methylation")==0){
			i++;
			output_meth = args[i];
			i++;
		}else{
			out("Error in the specified options.");
			help();
		}
		
		return i;
		
	}
	
	static int base_to_int(char c) {
		
		switch(c){
			case 'A': case 'a' : return 0;
			case 'C': case 'c' : return 1;
			case 'G': case 'g' : return 2;
			case 'T': case 't' : return 3;
		}
		
		return 4;
		
	}
	
	static char int_to_base(int i) {
		
		switch(i){
			case 0: return 'A';
			case 1: return 'C';
			case 2: return 'G';
			case 3: return 'T';
		}
		
		return 'A';
		
	}
	
	static boolean isNucleotide(char c){ return base_to_int(c)<4;}

	static char mutate(char c, String chr, long pos){
		
		if(!isNucleotide(c))
			return c;
		
		if(snp>0 && Math.random() < snp){
			
			int r = (((int)Math.random()*1000000)%3) + 1;
			c = int_to_base( (base_to_int(c) + r  )%4 );
			
		}
		
		//apply bisulfite conversion
		
		if(bs_fwd>0 && base_to_int(c)==base_to_int('C')){
			
			if(Math.random()<bs_fwd){
				
				c = 'T';
				
				try{
					if(output_meth.compareTo("")!=0) out_meth_annotation.write(chr+"\t"+pos+"\t"+pos+"\t+\t0\n");
				}catch(Exception e){}
				
			}else{
				
				try{
					if(output_meth.compareTo("")!=0) out_meth_annotation.write(chr+"\t"+pos+"\t"+pos+"\t+\t1\n");
				}catch(Exception e){}
				
			}
		
			
		}else if(bs_rev>0 && base_to_int(c)==base_to_int('G')){
			
			if(Math.random()<bs_rev){
				
				c = 'A';
				
				try{
					if(output_meth.compareTo("")!=0) out_meth_annotation.write(chr+"\t"+pos+"\t"+pos+"\t-\t0\n");
				}catch(Exception e){}
				
			}else{
				
				try{
					if(output_meth.compareTo("")!=0) out_meth_annotation.write(chr+"\t"+pos+"\t"+pos+"\t-\t1\n");
				}catch(Exception e){}
				
			}
			
		}
		
		
		return c;
		
	}
	
	//extract random indel length with exponential distribution
	static int randomIndelLength(){
		
		int l=1;
		
		while(Math.random()<extend)
			l++;
		
		return l;
		
	}
	
	static char randomNucleotide(){
		
		int r = ((int)(Math.random()*1000000)) % 4;
		return int_to_base(r);
		
	}
	
	public static String [] insertIndels(String seq,  String qual){
		
		String [] out = new String[2];
		out[0] = "";
		out[1] = "";
		int i=0;//position in original string		
		
		while(i<seq.length()){
		
			if(Math.random()<open){//open an indel here
				
				int l = randomIndelLength();
				
				if(Math.random()<0.5){//insertion in read
										
					for(int j=0;j<l;j++){
						out[0] += randomNucleotide();
						out[1] += qual.charAt(i);//copy current quality
					}
					
					//after the insertion, copy the character in the original read
					out[0] += seq.charAt(i);
					out[1] += qual.charAt(i);
					
					i++;
					
				}else{//deletion in read

					i+=l;//skip l bases in the original read
					
				}
				
			}else{
				out[0] += seq.charAt(i);
				out[1] += qual.charAt(i);			
				i++;
			}
			
		}
		
		return out;
		
	}
	
	public static String bsConvert(String seq){
		
		String out = "";
		
		for(int i=0;i<seq.length();i++){
			
			if(seq.charAt(i)=='C' || seq.charAt(i)=='c'){
				
				if(Math.random() < cr)
					out += 'T';
				else
					out += 'C';

			}else
				out += seq.charAt(i);
			
		}
		
			
		return out;
		
	}
	
	static void help(){
		out("Fastx Mutate Tools by Nicola Prezza: a set of tools to introduce simple mutations in fasta/fastq files\n");
		out("usage: 	java -jar FastxMutateTools.jar <snp|indel|bs> [options]");
		out("\n* snp:	introduce SNPs/bisulfite conversion in a fasta file. Options:");
		out("   --input = input fasta file. If not specified, standard input is used.");
		out("   --output = output fasta file. If not specified, standard output is used.");
		out("   --bs-fwd <freq> : convert a C into a T (bisulfite conversion on forward strand) with probability freq (in [0,1]). Default: 0.");
		out("   --bs-rev <freq> : convert a G into an A (bisulfite conversion on reverse strand) with probability freq (in [0,1]). Default: 0.");
		out("   --snp <freq> : introduce a SNP with probability freq (in [0,1]). If --bs-fwd or --bs-rev are specified, a nucleotide is firstly mutated with a snp, and then bisulfite-converted. Default: 0.");
		out("   --output-methylation <output-file> : output a file with the methylation values for each C/G. Format is: a line \"<contig> <pos> <pos> <strand> <beta>\" for each C and G in the input file. pos is the position of the C/G (repeated to have a well-formed bed file). strand is + or - if the conversion is C->T or G->A, respectively. beta is 0 (conversion introduced) or 1 (C/G not converted)");
		out("   --1-based : coordinates in the file specified with --output-methylation start from 1. Default: coordinates start from 0.");
		out("\n* indel: insert random indels in fastq file. Options:");
		out("   --input = input fastq file. If not specified, standard input is used.");
		out("   --output = output fastq file. If not specified, standard output is used.");
		out("   --open <freq> : probability to open an indel (in [0,1]). Default: 0.005.");
		out("   --extend <freq> : probability to extend an indel (in [0,1]). Default: 0.2.");
		out("\n* bs: insert random C->T bisulfite conversions in fastq file. Options:");
		out("   --input = input fastq file. If not specified, standard input is used.");
		out("   --output = output fastq file. If not specified, standard output is used.");
		out("   --cr <freq> : conversion rate (in [0,1]). Default: 0.8.");
		System.exit(0);
	}
	
	static void snpSimulator(){
	
		String chr="";//current chromosome name
		long pos = start_coord;
		
		try{
			
			while(in_fastx.ready()){
			
				char c = (char)in_fastx.read();
				
				if(c=='>'){
					
					String line = in_fastx.readLine();
					out_fastx.write(c+line+"\n");
					chr = line;
					pos = start_coord;
					
				}else{
				
					out_fastx.write(mutate(c,chr,pos));
					
				}
								
			}
			
			in_fastx.close();
			out_fastx.close();
			
			if(output_meth.compareTo("")!=0)
				out_meth_annotation.close();
			
		}catch(IOException e){
			out("Input/Output error:");
			e.printStackTrace();
		}
		
	}
	
	static void indelSimulator(){
		
		try{
			
			while (in_fastx.ready()) {
				
				out_fastx.write(in_fastx.readLine()+"\n");//name of read
								  
				String seq = in_fastx.readLine();
				in_fastx.readLine();//skip '+'
				String qual = in_fastx.readLine();
			  
				String [] seq_and_qual = insertIndels(seq,qual);
			  
				out_fastx.write(seq_and_qual[0]+"\n");
			  
				out_fastx.write("+\n");//'+'
			  
				out_fastx.write(seq_and_qual[1]+"\n");//quality
				  
			}
			
			in_fastx.close();
			out_fastx.close();
			
		}catch(Exception e){
			System.out.println("Error while reading file (maybe this is not a well-formatted fastq file?):");
			e.printStackTrace();
			System.exit(1);
		}
		
	}
	
	static void bsSimulator(){
		
		try{
			
			while (in_fastx.ready()) {
				
				out_fastx.write(in_fastx.readLine()+"\n");//name of read
				out_fastx.write(bsConvert(in_fastx.readLine())+"\n");//bs-converted sequence
				out_fastx.write(in_fastx.readLine()+"\n");//+
				out_fastx.write(in_fastx.readLine()+"\n");//quality
				  
			}
			
			in_fastx.close();
			out_fastx.close();
			
		}catch(Exception e){
			System.out.println("Error while reading file (maybe this is not a well-formatted fastq file?):");
			e.printStackTrace();
			System.exit(1);
		}
		
	}
	
	public static void main(String[] args) {
		
		if(args.length<1)//print help
			help();
				
		if(args[0].compareTo("snp")==0){
			mode=program_mode.snp_sim;
		}else if(args[0].compareTo("indel")==0){
			mode=program_mode.indel_sim;
		}else if(args[0].compareTo("bs")==0){
			mode=program_mode.bs_sim;
		}else{
			out("unrecognized program mode: " + args[0]);
			help();
		}
		
		int i = 1;
		
		while(i<args.length)
			i=parse(args,i);
		
		initBuffers();
		
		if(mode==program_mode.snp_sim)
			snpSimulator();
		
		if(mode==program_mode.indel_sim)
			indelSimulator();
		
		if(mode==program_mode.bs_sim)
			bsSimulator();
		
	}
	
}
