// problem-specific derived variables would be put here
//
// height
//
derive_lst.add("height",IndexType::TheCellType(),1,ca_derheight,the_same_box);
derive_lst.addComponent("height",desc_lst,State_Type,Density,NUM_STATE);

// primitive velocity
derive_lst.add("primu",IndexType::TheCellType(),1,ca_derprimu,the_same_box);
derive_lst.addComponent("primu",desc_lst,State_Type,Density,NUM_STATE);

derive_lst.add("primv",IndexType::TheCellType(),1,ca_derprimv,the_same_box);
derive_lst.addComponent("primv",desc_lst,State_Type,Density,NUM_STATE);

derive_lst.add("primw",IndexType::TheCellType(),1,ca_derprimw,the_same_box);
derive_lst.addComponent("primw",desc_lst,State_Type,Density,NUM_STATE);

// primitive density
derive_lst.add("primrho",IndexType::TheCellType(),1,ca_derprimrho,the_same_box);
derive_lst.addComponent("primrho",desc_lst,State_Type,Density,NUM_STATE);

// Lorentz factor
derive_lst.add("W",IndexType::TheCellType(),1,ca_derW,the_same_box);
derive_lst.addComponent("W",desc_lst,State_Type,Density,NUM_STATE);

// Internal energy
derive_lst.add("eint",IndexType::TheCellType(),1,ca_dereint,the_same_box);
derive_lst.addComponent("eint",desc_lst,State_Type,Density,NUM_STATE);

// Primitive momentum magnitude
derive_lst.add("primspeed",IndexType::TheCellType(),1,ca_derprimspeed,the_same_box);
derive_lst.addComponent("primspeed",desc_lst,State_Type,Density,NUM_STATE);

// Pressure gradient magnitute
derive_lst.add("pressgrad",IndexType::TheCellType(),1,ca_derpressgrad,the_same_box);
derive_lst.addComponent("pressgrad",desc_lst,State_Type,Density,NUM_STATE);
