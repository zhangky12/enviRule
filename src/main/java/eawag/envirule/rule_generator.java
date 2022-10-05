package eawag.envirule;

public class rule_generator {

    private boolean generalizeIgnoreHydrogen;
    private boolean includeFunctionalGroups;
    private String file;
    private int radius;

    public rule_generator(boolean generalizeIgnoreHydrogen, boolean includeFunctionalGroups, String file, int radius){
        this.generalizeIgnoreHydrogen = generalizeIgnoreHydrogen;
        this.includeFunctionalGroups = includeFunctionalGroups;
        this.file = file;
        this.radius = radius;
    }


}
