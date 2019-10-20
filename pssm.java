import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;

public class pssm {

	public static void main(String[] args) throws IOException {
		//1. read in the input sequence file + the motif training set
			//a. read in the sequence data.
			String input_file_name = args[0];
			//String input_file_name = "C:\\Users\\Owner\\Documents\\workspace\\PSSM\\src\\ecoK12-MG1655.fasta";
			FileReader input_seq_reader = new FileReader(input_file_name);
			BufferedReader input_seq_buff = new BufferedReader(input_seq_reader);
			StringBuilder input_seq = new StringBuilder("");
			String input_seq_line = input_seq_buff.readLine();
			while(input_seq_line != null) {
				if(input_seq_line.isEmpty())
					break;
				else if(input_seq_line.startsWith(">")) {
					input_seq_line = input_seq_buff.readLine();
					continue;
				}
				else {
					input_seq.append(input_seq_line);
					input_seq_line = input_seq_buff.readLine();
					continue;
				}
			}
			input_seq_buff.close();
			
				//TESTING
				//System.out.println("INPUT SEQUENCE SIZE: "+input_seq.length());
			//b. read in the motif training set
			String motif_file_name = args[1];
			//String motif_file_name = "C:\\Users\\Owner\\Documents\\workspace\\PSSM\\src\\FruR.fasta";
			FileReader motif_seq_reader = new FileReader(motif_file_name);
			BufferedReader motif_seq_buff = new BufferedReader(motif_seq_reader);
			StringBuilder motif_seq = new StringBuilder("");
			String motif_seq_line = motif_seq_buff.readLine();
			while(motif_seq_line != null) {
				if(motif_seq_line.isEmpty())
					break;
				else if(motif_seq_line.startsWith(">")) {
					motif_seq_line = motif_seq_buff.readLine();
					continue;
				}
				else {
					motif_seq.append(motif_seq_line);
					motif_seq.append(":");
					motif_seq_line = motif_seq_buff.readLine();
					continue;
				}
			}
			motif_seq_buff.close();
			String[] motif_set = motif_seq.toString().split(":");
			
				//TESTING
				//for(String s: motif_set) {System.out.println(s);}
				//System.out.println(motif_set.length);
			
		//2. create the pssm
			//a. count the nucleotide frequencies at each site, score the positions.
			double[][] pssm = new double[motif_set[0].length()][4];
			int a = 0;
			int t = 0;
			int c = 0;
			int g = 0;
			for(int i = 0; i < motif_set[0].length();i++) {
				for(int j = 0; j < motif_set.length;j++) {
					switch(motif_set[j].toUpperCase().charAt(i)) {
					case 'A':
						a++;
						break;
					case 'T':
						t++;
						break;
					case 'C':
						c++;
						break;
					case 'G':
						g++;
						break;
					}
				}
				//* A higher score means that the nucleotide is the most frequent.
				pssm[i][0] = Math.log(((a + 0.25)/motif_set.length)/0.25);
				pssm[i][1] = Math.log(((t + 0.25)/motif_set.length)/0.25);
				pssm[i][2] = Math.log(((c + 0.25)/motif_set.length)/0.25);	
				pssm[i][3] = Math.log(((g + 0.25)/motif_set.length)/0.25);	
				
				a = 0;
				t = 0;
				c = 0;
				g = 0;
			}
			
			//b. suggest a reasonable cutoff for future studies
			
		
			System.out.println("MOTIF TRAINING SET PSSM:");
				for(double[] d: pssm) {
					System.out.println("\t"+Arrays.toString(d));
				}
			System.out.println();
			System.out.println("SCORES FOR MOTIF TRAINING DATA:");
				double indiv_score = 0.0;
				for(int i = 0;i <  motif_set.length;i++) {
					for(int j = 0;j < motif_set[0].length();j++) {
						switch(motif_set[i].toUpperCase().charAt(j)) {
						case 'A':
							indiv_score += pssm[i][0];
							break;
						case 'T':
							indiv_score += pssm[i][1];
							break;
						case 'C':
							indiv_score += pssm[i][2];
							break;
						case 'G':
							indiv_score += pssm[i][3];
							break;
						}
					}
					System.out.println("\t"+motif_set[i]+" "+indiv_score);
					indiv_score = 0.0;
				}
			System.out.println();
			
		//3. scan the input sequence  for motifs, their positions, and their strand (+,-)
				//a. implement a sliding-window approach. initialize cutoff point
				double cutoff = Double.parseDouble(args[2]);
				//double cutoff = 0.0;
				int window_start = 0;
				double score = 0.0;
				double plus_score = 0.0;
				double minus_score = 0.0;
				StringBuilder window = new StringBuilder("");
				String strand = "+";
				System.out.println("MOTIFS IN INPUT SEQUENCE");
				while(window_start + motif_set[0].length() - 1 < input_seq.length()) {
					window.append(input_seq.substring(window_start, window_start + motif_set[0].length()).toUpperCase());
					for(int i = 0;i < window.length();i++) {
						switch(window.charAt(i)) {
						case 'A':
							plus_score += pssm[i][0];
							break;
						case 'T':
							plus_score += pssm[i][1];
							break;
						case 'C':
							plus_score += pssm[i][2];
							break;
						case 'G':
							plus_score += pssm[i][3];
							break;
						}
						
						switch(window.reverse().charAt(i)) {
						case 'A':
							minus_score += pssm[i][0];
							break;
						case 'T':
							minus_score += pssm[i][1];
							break;
						case 'C':
							minus_score += pssm[i][2];
							break;
						case 'G':
							minus_score += pssm[i][3];
							break;
						}
					}
					
					score = plus_score;
					if(plus_score < minus_score) {
						strand = "-";
						score = minus_score;
						window = window.reverse();
					}
					
					if(score > cutoff) {
						System.out.println("\t"+window+" "+score+" "+strand+" "+window_start+" "+(window_start+motif_set[0].length()-1));
					}
					window.setLength(0);
					window_start++;
					score = 0.0;
					plus_score = 0.0;
					minus_score = 0.0;	
					strand = "+";
				}

				
			}
			
}
