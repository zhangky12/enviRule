package eawag.envirule;

import eawag.envirule.servelet.rule_generator_servelet;
import eawag.envirule.servelet.rxn_adder_servelet;
import eawag.envirule.servelet.rxn_clusterer_servelet;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.ObjectOutputStream;
import java.net.ServerSocket;
import java.net.Socket;
import java.util.Set;

public class enviRule_server {

    public enviRule_server(int port) {

        //initialize socket and input stream
        Socket socket;
        ServerSocket server;

        // starts server and waits for a connection
        try{
            server = new ServerSocket(port);
            server.setReuseAddress(true);

            System.out.println("Welcome to enviRule developed by EAWAG!");
            System.out.println("Waiting for a client ...");

            while (true){
                socket = server.accept();
                System.out.println("New client connected "
                        + socket.getInetAddress()
                        .getHostAddress());
                enviRule_server.ClientHandler clientSock = new enviRule_server.ClientHandler(socket);
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
                            res = rule_generator_servelet.test(line.split(" "));
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
                            rxn_clusterer_servelet.cluster(line.split(" "));
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
                            res = rxn_adder_servelet.adder(line.split(" "));
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
        enviRule_server server = new enviRule_server(5000);
    }
}
