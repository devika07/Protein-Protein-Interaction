import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.Map.Entry;
import java.util.Queue;
import org.graphstream.*;
import org.graphstream.algorithm.*;
import org.graphstream.graph.Graph;
import org.graphstream.graph.Node;
import org.graphstream.graph.implementations.MultiGraph;
public class Protein_analysis {
	static LinkedHashMap<String,ArrayList<String>> hk;
	static LinkedHashMap<String,Integer> ht;
	static LinkedHashMap<String,Integer> hd;
	static int diameter;     //diameter of the graph
	static int degree;       //degree of each node
	static int radius;       //radius of the graph	 
	@SuppressWarnings("unused")
	public static void main(String[] args) throws IOException
	{
		String filename="fly_net.txt";
		BufferedReader br = new BufferedReader(new FileReader(filename));
		String line;
		hk=new LinkedHashMap<String, ArrayList<String>>();
		ht=new LinkedHashMap<String, Integer>();
		hd=new LinkedHashMap<String,Integer>();
		while ((line = br.readLine()) != null)
		{
			String[] tokens=line.split("\\:|\\ ");
			String P1=tokens[1];   //Protein 1
			String P2=tokens[3];   //Protein 2
			if(!(P1.equals(P2))){
			if(!(hk.containsKey(P2))){
				ArrayList<String> m=new ArrayList<String>();
				m.add(P1);
				hk.put(P2,m);}
			else{
				ArrayList<String> o=new ArrayList<String>();
				o=hk.get(P2);
				o.add(P1);
				hk.put(P2, o);}
			if(!(hk.containsKey(P1))){
				ArrayList<String> a=new ArrayList<String>();
				a.add(P2);
				hk.put(P1,a);}
			else{ 
		     ArrayList<String> b=new ArrayList<String>();
		     b=hk.get(P1);
		     b.add(P2);
		     hk.put(P1, b);}
			}
		}
		HashMap<String,Integer> ht=new LinkedHashMap<String,Integer>();
		for (Entry<String, ArrayList<String>> entry : hk.entrySet()) {
		    String key = entry.getKey();
		    ArrayList<String> value = entry.getValue();
		        int x=breadthFirstSearch(key);
		        ht.put(key, x);}
		diameter=0;
		for(Entry<String,Integer>entry1:ht.entrySet()){
			int y=entry1.getValue();
			if (y>diameter)
				diameter=y;}   //diameter of graph
		radius=1000; 
		for(Entry<String,Integer>entry2:ht.entrySet()){
			int y=entry2.getValue();
			if (y<radius)
				radius=y;}      //radius of graph
		int a=degree_max();
		double z=clustering_coefficient();
System.out.println("The diameter is : "+ diameter);
System.out.println("The maximum degree Centrality is : "+ a	);
System.out.println("The average clustering coefficient is : "+z);
System.out.println("The radius is :"+radius);
//Adding nodes to graph
@SuppressWarnings("deprecation")
Graph graph=new MultiGraph("Graph1");
ArrayList<Node> array=new ArrayList<Node>();
for (Entry<String, ArrayList<String>> entry : hk.entrySet()) {
	String key = entry.getKey();
	Node n=graph.addNode(key);
	array.add(n);}
for (Entry<String, ArrayList<String>> entry : hk.entrySet()) {
	String key4 = entry.getKey();
	ArrayList<String> value = entry.getValue();
	for(String nodestring:value){
		String g=key4+"_"+nodestring;
		String h=nodestring+"_"+key4;
		if(!(graph.hasLabel(g)||graph.hasLabel(h)))
		graph.addEdge(key4+"_"+nodestring,key4,nodestring);}
}
//Betweenness Centrality
BetweennessCentrality bcb = new BetweennessCentrality();
bcb.init(graph);
bcb.compute();
double max=0;
for(Node s:array)
{
   double k=s.getAttribute("Cb");
	if(k>max)
		max=k;
}
int b=hk.size();
max=max/((b-1)*(b-2));     //Betweenness Centrality
System.out.println("The maximum of betweenness centrality is :"+max);
graph.display();           //displaying the graph
}
	//BFS to find the eccentricity of each node
		public static int breadthFirstSearch(String key){
		HashMap<String,Integer> hl=new LinkedHashMap<String,Integer>();
			Queue<String> q = new LinkedList<String>();
			int e=0;
			hl.put(key,0);
			q.add(key);
			while(!q.isEmpty()){
				String u = q.poll();
				if(hk.containsKey(u)){
				for (String entry : hk.get(u)){
					if(!(hl.containsKey(entry))){
						int w= hl.get(u)+1;
						hl.put(entry,w);
						e = w;
						q.add(entry);}
					else if(entry.equals(u)){
						int w=hl.get(u)+1;
						hl.put(entry,w);
						e=w;}
				}}}
			return e;
		}
//Finding the degree of each node
public static int degree_max(){
	int max=0;
	for(Entry<String,ArrayList<String>>entry1:hk.entrySet()){
		ArrayList<String> value = entry1.getValue();
		hd.put(entry1.getKey(), value.size());
		 if(value.size()>max)
			 max=value.size();
	}
	return max;
}
//Calculating Local Clustering Coefficient for each node and then finding average
public static double clustering_coefficient()
{
	double avg=0;int k;int t;int f=0;
for(Entry<String,ArrayList<String>>entry2:hk.entrySet()){
		t=0;
		f++;
		ArrayList<String> c=new ArrayList<String>();
	    c=entry2.getValue();
	k=hd.get(entry2.getKey());
	for(int i=0;i<c.size();i++)
	{
		String a=c.get(i);
		for(int j=i+1;j<c.size();j++)
		{
			String b=c.get(j);
			if(hk.containsKey(a)){
				if(hk.get(a).contains(b))
					t++;}
		}}
	if(k!=1)
	avg=avg+(2.0*t/(k*k-k));
}
	return avg/f;
}
}