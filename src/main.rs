use bevy::prelude::*;
use bevy_rapier2d::prelude::*;
use components::camera_movement;
use components::check_sv;
use components::confine_camera_movement;
use components::confine_player_movement;
use components::player_movement;
use components::set_jumping_false_if_touching_floor;
use components::spawn_camera;
use components::spawn_landscape;
use components::spawn_player;

pub mod grid;
pub mod components;


fn main() {

    App::new()
        .add_plugins((DefaultPlugins, RapierPhysicsPlugin::<NoUserData>::default()))
        .add_systems(Startup, (spawn_camera, spawn_landscape, spawn_player.after(spawn_landscape)))
        .add_systems(Update, (player_movement,
                              set_jumping_false_if_touching_floor,
                              camera_movement.after(player_movement),
                              check_sv,
                              //confine_player_movement.after(player_movement),
                              confine_camera_movement))
        .run();

}
