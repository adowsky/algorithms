package algorithm.graphs;


import algorithm.graphs.heuristics.MarkovClustering;
import algorithm.graphs.representation.WeightenedIncidenceMatrix;

import java.io.*;
import java.util.LinkedList;

/**
 * Reads DNA sequences from file and creates graph with subseqences.
 */
public class DNAPatterns {
    private WeightenedIncidenceMatrix<DNANode> graph;
    private int hamm;
    private int instancesNumber;
    private int nodeLength;
    private boolean devel;
    public DNAPatterns(String fasta,String qual, int length, int hamm){
        this.hamm = hamm;
        this.nodeLength = length;
        try(BufferedReader inFasta = new BufferedReader(new FileReader(new File(fasta)));
                InputStreamReader inQual = new InputStreamReader(new BufferedInputStream(new FileInputStream(qual)))){
            String line;
            int seqNo = 0;
            int pos = 0;
            DNANode[] nodes;
            DNANode[][] nodesArrs = new DNANode[10][];
            while((line = inFasta.readLine()) != null){
                if(line.charAt(0) == '>'){//Reading sequence identification
                    String[] chopped = line.split(" ");
                    ++seqNo;
                    int seqLen = Integer.valueOf(chopped[1].split("=")[1]);
                    line = inFasta.readLine();
                    line = line.concat(inFasta.readLine());
                    DNANode[] current = new DNANode[seqLen-length+1];
                    for(int i=0; i<(seqLen- length +1);i++){
                        String part = line.substring(i,i+length);
                        current[i] = new DNANode(length,part,seqNo,i);
                    }
                    nodesArrs[seqNo-1] = current;
                }
            }
            instancesNumber = seqNo;
            int len = 0;
            for(DNANode[] e :nodesArrs){
                if(e!=null)
                    len += e.length;
            }
            nodes = new DNANode[len];
            int iter =0;
            for(DNANode[] e :nodesArrs){
                if(e!= null)
                for(DNANode n : e){
                    nodes[iter++] = n;
                }
            }
            graph = new WeightenedIncidenceMatrix<>(nodes);

        }catch(IOException ex){
            ex.printStackTrace();
        }
    }
    public void findConnections(){
        if(devel)
            System.out.print("Making connections between vertices... ");
        for(int i=0;i<graph.size();++i) {
            DNANode node = graph.getVertex(i);
            for (int j = 0; j < graph.size(); ++j)
                if (i != j && node.getInstNo()!=graph.getVertex(j).getInstNo() && node.isSimiliarTo(graph.getVertex(j)) ){
                    int k = node.hamming(graph.getVertex(j));
                    if(k <= hamm){
                        k = 10*(node.length - k);
                    }else {
                        k = node.length - k;
                    }
                    graph.setConnection(i,j,k);
                    graph.setConnection(j,i,k);
                }
        }

        if(devel)
            System.out.println("[OK]");
    }

    private void findResult(int[][] cls, int size){
        LinkedList<DNANode>[] tmp= new LinkedList[instancesNumber];
        StringBuilder[] values = new StringBuilder[instancesNumber];
        LinkedList<DNANode>[] result= new LinkedList[instancesNumber];
        StringBuilder resultVal = null;
        for(int i=0;i<instancesNumber;++i) {
            tmp[i] = new LinkedList<>();
            values[i] = new StringBuilder();
        }
            for(int[] clique : cls){
                if(clique == null || clique.length < 2)
                    continue;
                for(int index : clique){
                    DNANode vertex = graph.getVertex(index);
                    if(tmp[vertex.instNo -1 ].size() > 0) {
                        DNANode last = tmp[vertex.getInstNo() - 1].getLast();
                        if (vertex.getPos() - last.getPos() < vertex.getLength() && vertex.getPos() - last.getPos() >0) {
                            tmp[vertex.instNo - 1].addLast(vertex);
                            String str = vertex.getSeq().substring(vertex.getLength() -(vertex.getPos() - last.getPos()));
                            values[vertex.instNo - 1].append(str);
                        } else {
                            if (resultVal == null) {
                                resultVal = values[vertex.instNo - 1];
                            }
                            values[vertex.instNo - 1] = new StringBuilder();
                            tmp[vertex.instNo - 1] = new LinkedList<>();
                            tmp[vertex.instNo - 1].addLast(vertex);
                        }
                    }else
                        tmp[vertex.instNo -1 ].addLast(vertex);
                }
        }
        System.out.println("Best Theme: "+resultVal.toString());
    }

    public void setDeveloper(boolean t){
        devel = t;
    }
    public static void main(String[] args){
        DNAPatterns pt = new DNAPatterns("s1.fast", "s1.qual",4,1);
        pt.setDeveloper(true);
        pt.findConnections();
        MarkovClustering mc = new MarkovClustering(pt.graph.getMatrix(), 2,2);
        mc.setDevel(true);
        pt.findResult(mc.execute(),mc.getCliquesSize());

    }
}
class DNANode{
    int length;
    String seq;
    int instNo;
    int pos;

    public DNANode(int l, String seq, int instNo, int pos ){
        length = l;
        this.pos = pos;
        if(seq.length() != length)
            System.out.println("Length Error!");
        this.seq = seq;
        this.instNo = instNo;
    }

    public int getPos() {
        return pos;
    }

    public int getLength() {
        return length;
    }

    public String getSeq() {
        return seq;
    }

    public int getInstNo() {
        return instNo;
    }

    @Override
    public boolean equals(Object o){
        //TODO need to be improved
       return ((DNANode)o).getSeq().equals(seq);
    }
    public boolean isSimiliarTo(DNANode o){
        int k=0;
        for(int i=0;i<length;i++){
            if(seq.charAt(i) != o.seq.charAt(i))
                k++;
        }
        if(k>2)
            return false;
        else
            return true;
    }
    public int hamming(DNANode o){
        int r = 0;
        for(int i=0;i<length;++i){
            if(seq.charAt(i) != o.getSeq().charAt(i))
                r++;
        }
        return r;
    }
}
class ResultStruct{
    private DNANode[] sq;
    private int[][] conns;
    private int size;
    private int[] vertices;

    public int[] getVertices() {
        return vertices;
    }

    public void setVertices(int[] vertices) {
        this.vertices = vertices;
    }

    public DNANode[] getSq() {
        return sq;
    }

    public ResultStruct setSq(DNANode[] sq) {
        this.sq = sq;
        return this;
    }

    public int[][] getConns() {
        return conns;
    }

    public ResultStruct setConns(int[][] conns) {
        this.conns = conns;
        return this;
    }

    public int getSize() {
        return size;
    }

    public ResultStruct setSize(int size) {
        this.size = size;
        return this;
    }
}
class LongestPathFinder{
    private int[][] graph;
    public  LongestPathFinder(int[][] graph){
        this.graph = graph;
    }
    public int[] getLongestPath(){
        int[] result = new int[graph.length];
        int[] tmp = new int[graph.length];
        int iterator = 0;
        boolean wasNext = true;
        int bestValue = 0;
        int value = 0;
        //do we need to check start from each vertex?
        for(int k =0; k<graph.length;++k) {
            tmp[iterator++] = k;
            do {
                int v = tmp[iterator - 1];
                int start = 0;
                if (!wasNext) {
                    start = tmp[iterator];
                    tmp[iterator] = -1;
                }
                wasNext = false;
                for (int i = start; i < graph.length; ++i) {
                    if (graph[v][i] > 0) {
                        value += graph[v][i];
                        tmp[iterator++] = i;
                        wasNext = true;
                        break;
                    }
                }
                if(!wasNext){
                    if(value> bestValue){
                        for(int i=0;i<graph.length;++i){
                            result[i] = tmp[i];
                            bestValue = value;
                        }
                    }
                    iterator--;
                }

            } while (iterator>0);
        }
        return result;
    }

}