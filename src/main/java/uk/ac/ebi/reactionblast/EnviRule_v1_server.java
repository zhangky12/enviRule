package uk.ac.ebi.reactionblast;
import java.net.*;
import java.io.*;
import java.text.ParseException;
import java.util.Set;

public class EnviRule_v1_server {

    //initialize socket and input stream
    private Socket          socket   = null;
    private ServerSocket    server   = null;

    // constructor with port
    public EnviRule_v1_server(int port) {
        // starts server and waits for a connection
        try{
            server = new ServerSocket(port);
            server.setReuseAddress(true);

            System.out.println("Server started");

            System.out.println("Waiting for a client ...");

            while (true){
                socket = server.accept();
                System.out.println("New client connected "
                        + socket.getInetAddress()
                        .getHostAddress());
                ClientHandler clientSock = new ClientHandler(socket);
                System.out.println("New thread created");
                new Thread(clientSock).start();
            }
        }
        catch(IOException i){
            System.out.println(i);
        }
    }

    private static class ClientHandler implements Runnable {
        private final Socket clientSocket;

        public ClientHandler(Socket socket){
            this.clientSocket = socket;
        }

        @Override
        public void run() {

            ObjectOutputStream out = null;
            BufferedReader in = null;
            try{
                out = new ObjectOutputStream(clientSocket.getOutputStream());
                in = new BufferedReader(new InputStreamReader(clientSocket.getInputStream()));

                String line;
                while ((line = in.readLine()) != null){

                    if(line.length() < 9){
                        System.out.println("Invalid command");
                        out.writeObject(null);
                        continue;
                    }

                    if(line.substring(0, 9).compareTo("autorule ")==0){
                        System.out.println("Generating rules...");
                        Set<String> res = null;
                        try{
                            res = EnviRule_v1.test(line.split(" "));
                        }catch (Exception e){
                            System.out.println(e.getMessage());
                        }

                        if(res == null) {
                            out.writeObject(null);
                            continue;
                        }
                        out.writeObject(res);

                    }else if(line.substring(0, 9).compareTo("rxnclust ")==0){
                        System.out.println("Clustering reactions...");
                        boolean success = false;
                        try{
                            RxnCluster_v1.cluster(line.split(" "));
//                            String[] input = line.split(" ");
//                            Clustering_rxn_class.clustering(input[1], input[2]);
                            success = true;
                        }catch (Exception e){
                            System.out.println(e.getMessage());
                        }

                        if(!success){
                            out.writeObject(null);
                            continue;
                        }
                        out.writeObject("Clustering done");
                    }else if(line.substring(0, 9).compareTo("rxnadder ")==0){
                        System.out.println("Adding new reactions...");
                        Set<String> res = null;
                        boolean success = false;
                        try{
                            res = RxnAdder_v1.adder(line.split(" "));
                            success = true;
                        }catch (Exception e){
                            System.out.println(e.getMessage());
                        }

                        if(!success){
                            out.writeObject(null);
                            continue;
                        }
                        out.writeObject(res);
                    }else{
                        System.out.println("Invalid command");
                        out.writeObject(null);
                    }
                }
            } catch (Exception e) {
                e.printStackTrace();
            } finally {
                System.out.println("Closing connection...");
                try{
                    if (out != null) out.close();
                    if (in != null) {
                        in.close();
                        clientSocket.close();
                    }
                }catch (IOException e){
                    e.printStackTrace();
                }
            }
            System.out.println("Thread closed");
        }
    }

    public static void main(String args[]){
        EnviRule_v1_server server = new EnviRule_v1_server(5000);
    }


}
