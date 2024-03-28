use bevy::prelude::*;
use components::camera_movement;
use components::confine_camera_movement;
use components::confine_player_movement;
use components::player_movement;
use components::spawn_camera;
use components::spawn_landscape;
use components::spawn_player;

pub mod grid;
pub mod components;
use std::vec;

use crate::grid::Direction;
use crate::grid::Grid;
use crate::grid::Player;
use crate::grid::get_input;

fn main() {

    App::new()
        .add_plugins(DefaultPlugins)
        .add_startup_system(spawn_camera)
        .add_startup_system(spawn_landscape)
        .add_startup_system(spawn_player)
        .add_system(player_movement)
        .add_system(camera_movement)
        .add_system(confine_player_movement)
        .add_system(confine_camera_movement)
        .run();

}
