def write_xyz(filename, positions, comment=""):
    with open(filename, "w") as f:
        f.write(f"{len(positions)} \n")
        f.write(f"{comment} \n")
        for i in positions:
            f.write(f"C    {i[0]:<20}  {i[1]:<20}  {i[2]:<20} \n")
            
            
            
def add_xyz_frame(filename, positions, comment=""):
    with open(filename, "a") as f:
        f.write(f"{len(positions)} \n")
        f.write(f"{comment} \n")
        for i in positions:
            f.write(f"C    {i[0]:<20}  {i[1]:<20}  {i[2]:<20} \n")       